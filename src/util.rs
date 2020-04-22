use std::cmp::Ordering;
use std::collections::{HashMap, HashSet};
use std::convert::TryInto;
use std::f64;
use std::fs::*;
use std::io::prelude::*;
use std::io::{self, BufReader, BufWriter};
use std::time::Instant;
//use std::num::pow::pow;
//use std::collections::BinaryHeap;

use math::round;
use binary_heap_plus::*;
use byteorder::{ByteOrder, LittleEndian};
use flate2::read::GzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;
use itertools::Itertools;
use nalgebra as na;
use ndarray::prelude::*;
use num_traits::cast::ToPrimitive;
use ordered_float::*;
use petgraph as pg;
use petgraph::unionfind::UnionFind;
use petgraph::visit::EdgeRef;
use rand::distributions::{Distribution, Uniform};
use rand_core::SeedableRng;
use rand_pcg::Pcg64;
use refinery::Partition;
use serde::Deserialize;
//use rand::thread_rng;
//use rgsl::statistics::correlation;

use crate::salmon_types::{AlevinMetaData, EdgeInfo, EqClassExperiment, BFHEqClassExperiment, FileList, MetaInfo, TxpRecord};

// General functions to r/w files
// files to be handled
// quant.sf
// metainfo.json
// eq_classes.txt
// bootstraps.gz
// ambig_info.tsv

pub fn group_writer(
    gfile: &mut File,
    groups: &HashMap<usize, Vec<usize>>,
) -> Result<bool, io::Error> {
    //let mut buffer = File::create("groups.txt")?;
    //let mut buffer = File::create("foo.txt").unwrap();
    //let mut file = GzEncoder::new(file_handle, Compression::default());
    for (group_id, group) in groups {
        let strings: Vec<String> = group.iter().map(|n| n.to_string()).collect();
        writeln!(gfile, "{},{}", group_id.to_string(), strings.join(","))?;
    }
    Ok(true)
}

pub fn write_quants_from_components(
    components: &[Vec<pg::graph::NodeIndex>],
    file_list: &FileList,
    gibbs_mat: &Array2<f64>,
    meta_info: &MetaInfo,
    rec: &[TxpRecord],
) -> Result<bool, io::Error> {
    let aux_dir = file_list.eq_file.parent().unwrap();
    create_dir_all(aux_dir)?;
    let bootstraps_dir = file_list.bootstrap_file.parent().unwrap();
    create_dir_all(bootstraps_dir)?;

    // write bootstrap file
    let bt_file_handle = File::create(file_list.bootstrap_file.clone())?;
    let mut file = GzEncoder::new(bt_file_handle, Compression::default());

    // write quantfile
    let quant_file =
        File::create(file_list.quant_file.clone()).expect("Could not write quant file");
    let mut quant_buf_writer = BufWriter::new(quant_file);

    // write meta_info.jeson file
    let mut json_file = File::create(file_list.mi_file.clone()).expect("Could not write json file");

    // write names.tsv.gz file
    let name_file_handle = File::create(file_list.names_tsv_file.clone())?;
    let mut name_file = GzEncoder::new(name_file_handle, Compression::default());

    // write cluster file
    let mut cluster_file =
        File::create(file_list.cluster_file.clone()).expect("Could not write cluster file");

    let shape = gibbs_mat.shape();
    let new_height = components.len();
    let new_width = shape[1];

    // new structure for gibbs
    let mut gibbs_new_mat = Array2::<f64>::zeros((new_height, new_width));

    // new structure for quant
    // let rec = parse_quant(&quant_file, &meta_info);
    let mut quant_new_rec: Vec<TxpRecord> = Vec::new();
    quant_new_rec.resize_with(new_height, Default::default);

    // meta file
    // println!("{:?}, {}, {}", file_list.mi_file, new_height, shape[0]);
    let mut meta_info_new = meta_info.clone();
    meta_info_new.num_valid_targets = new_height.try_into().unwrap();
    let meta_info_data = serde_json::to_string_pretty(&meta_info_new).unwrap();
    // println!("{:?}", meta_info_data);
    json_file.write_all(meta_info_data.as_bytes())?;

    // names.tsv.gz
    let mut name_vec: Vec<String> = vec!["".to_string(); new_height];

    //let mut new_map = vec![0 as usize; new_height];
    //let mut ind = 0 as usize;

    for (i, comp) in components.iter().enumerate() {
        if comp.len() == 1 {
            // write it as it is
            let mut s = gibbs_new_mat.slice_mut(s![i, ..]);
            let component_index = comp[0].index();
            let to_initiate = gibbs_mat.index_axis(Axis(0), component_index).to_owned();
            s += &to_initiate;

            if let Some(txp_rec_elem) = quant_new_rec.get_mut(i) {
                let s = rec[component_index].clone();
                txp_rec_elem.assign(&s);
            }

            name_vec[i] = rec[component_index].Name.clone();
        } else {
            // collpase groups

            // the slice that will hold the gibbs samples for the group
            let mut s = gibbs_new_mat.slice_mut(s![i, ..]);

            // the sum of TPM for elements in the group (for normalizing lengths)
            let group_tpm = comp
                .iter()
                .map(|g| rec[g.index()].TPM)
                .fold(0_f64, |sum, e| sum + e as f64);
            let nonzero = group_tpm > 0_f64;

            // rewrite quants
            let new_name = format!("NewTr{}", i);
            name_vec[i] = new_name.clone();

            let mut member_names: Vec<String> = Vec::with_capacity(comp.len());

            if let Some(txp_rec_elem) = quant_new_rec.get_mut(i) {
                let mut weighted_len;
                let mut weighted_eff_len;

                let mut component_index = comp[0].index();

                s += &gibbs_mat.index_axis(Axis(0), component_index).to_owned();

                let initial_elem = rec[component_index].clone();
                member_names.push(initial_elem.Name.clone());

                txp_rec_elem.assign(&initial_elem);
                // the name is the new name
                txp_rec_elem.Name = new_name.clone();
                // the starting effective length is the *weighted* effective length
                if nonzero {
                    let expression_ratio = initial_elem.TPM as f64 / group_tpm;
                    weighted_eff_len = (initial_elem.EffectiveLength as f64) * expression_ratio;
                    weighted_len = (initial_elem.Length as f64) * expression_ratio;
                } else {
                    weighted_eff_len = initial_elem.EffectiveLength as f64;
                    weighted_len = initial_elem.Length as f64;
                }

                for g in comp.iter().skip(1) {
                    component_index = g.index();
                    s += &gibbs_mat.index_axis(Axis(0), component_index).to_owned();

                    let old_elem = rec[component_index].clone();
                    let expression_ratio = old_elem.TPM as f64 / group_tpm;
                    txp_rec_elem.TPM += old_elem.TPM;
                    txp_rec_elem.NumReads += old_elem.NumReads;
                    if nonzero {
                        weighted_len += (old_elem.Length as f64) * expression_ratio;
                        weighted_eff_len += (old_elem.EffectiveLength as f64) * expression_ratio;
                    } else {
                        weighted_eff_len = if old_elem.EffectiveLength as f64 > weighted_eff_len {
                            old_elem.EffectiveLength as f64
                        } else {
                            weighted_eff_len
                        };

                        weighted_len = if old_elem.Length as f64 > weighted_len {
                            old_elem.Length as f64
                        } else {
                            weighted_len
                        };
                    }
                    member_names.push(old_elem.Name);
                }

                // set the length by rounding the weighted length
                txp_rec_elem.Length = weighted_len.round() as u32;
                txp_rec_elem.EffectiveLength = weighted_eff_len as f32;
            }

            writeln!(cluster_file, "{},{}", new_name, member_names.join(","))?;
        }
    }

    // write bootstrap file
    gibbs_new_mat = gibbs_new_mat.reversed_axes();
    //println!("{:?}, {}", gibbs_new_mat.shape(),gibbs_new_mat.len());
    assert!([new_width, new_height] == gibbs_new_mat.shape());

    for i in 0..(new_width) {
        let mut data: Vec<u8> = vec![0; (new_height * 8) as usize];
        //br.read_exact(&mut data[..]).expect("could not read from inferential replicate file");
        //let mut floats: na::DVector<f64> = na::DVector::<f64>::from_element(nt, 0.0);// Vec<f64> = vec![0.0_f64; nt];
        LittleEndian::write_f64_into(&gibbs_new_mat.slice(s![i, ..]).to_vec(), &mut data);
        file.write_all(&data)?;
    }

    // write quant.sf file
    // println!("{}",msg);
    quant_buf_writer
        .write_all(b"Name\tLength\tEffectiveLength\tTPM\tNumReads\n")
        .expect("couldn't write to quant.sf");

    for elem in &quant_new_rec {
        quant_buf_writer
            .write_all(
                &format!(
                    "{}\t{}\t{}\t{}\t{}\n",
                    elem.Name, elem.Length, elem.EffectiveLength, elem.TPM, elem.NumReads
                )
                .as_bytes(),
            )
            .expect("couldn't write to quant.sf");
    }

    name_file.write_all(name_vec.join("\t").as_bytes())?;
    Ok(true)
}

#[allow(dead_code)]
pub fn parse_quant(p: &std::path::Path, mi: &MetaInfo) -> Result<Vec<TxpRecord>, io::Error> {
    let file = File::open(p);
    let mut rdr = csv::ReaderBuilder::new()
        .has_headers(true)
        .delimiter(b'\t')
        .from_reader(file.unwrap());

    let mut rs = Vec::<TxpRecord>::new();
    for result in rdr.deserialize() {
        // Notice that we need to provide a type hint for automatic
        // deserialization.
        let record: TxpRecord = result?;
        rs.push(record)
        //println!("{:?}", record);
    }
    assert!(
        rs.len() == mi.num_valid_targets as usize,
        "obsered wrong number of targets in quantification file!"
    );
    Ok(rs)
}

pub fn parse_json(p: &std::path::Path) -> Result<MetaInfo, io::Error> {
    let file = File::open(p);
    let reader = BufReader::new(file.unwrap());

    let jd: MetaInfo = serde_json::from_reader(reader).unwrap();
    println!("# targets : {}", jd.num_valid_targets);
    println!("did serialize eq classes : {}", jd.serialized_eq_classes);
    println!("# boot : {}", jd.num_bootstraps);
    Ok(jd)
}

#[allow(dead_code)]
pub fn read_gibbs(f: &std::path::Path, mi: &MetaInfo, gibbs_mat: &mut na::DMatrix<f64>) {
    let file = File::open(f);
    let s = GzDecoder::new(file.unwrap());
    let mut br = io::BufReader::new(s);

    let nt = mi.num_valid_targets as usize;
    for i in 0..(mi.num_bootstraps) {
        let mut data: Vec<u8> = vec![0; (nt * 8) as usize];
        br.read_exact(&mut data[..])
            .expect("could not read from inferential replicate file");
        let mut floats: na::DVector<f64> = na::DVector::<f64>::from_element(nt, 0.0); // Vec<f64> = vec![0.0_f64; nt];
        LittleEndian::read_f64_into(&data, floats.as_mut_slice());
        let mut col = gibbs_mat.column_mut(i as usize);
        col.copy_from(&(&col + floats));
    }

    //floats.extend_from_slice(unsafe{ std::mem::transmute::<&[u8], &[f64]>(&data[0..(nt*8)]) });
    //println!("{:?}", size);
}

pub fn read_gibbs_array(f: &std::path::Path, mi: &MetaInfo, gibbs_mat: &mut Array2<f64>) {
    let file = File::open(f);
    let s = GzDecoder::new(file.unwrap());
    let mut br = io::BufReader::new(s);

    let nt = mi.num_valid_targets as usize;
    for i in 0..(mi.num_bootstraps) {
        let mut data: Vec<u8> = vec![0; (nt * 8) as usize];
        br.read_exact(&mut data[..])
            .expect("could not read from inferential replicate file");
        let mut floats = vec![0.0; nt]; // Vec<f64> = vec![0.0_f64; nt];
        LittleEndian::read_f64_into(&data, &mut floats);
        gibbs_mat
            .slice_mut(s![.., i as usize])
            .assign(&Array::from(floats))
        //let mut col = gibbs_mat.column_mut(i as usize);
        //col.assign(&float);
    }

    //println!("Value at location 0,20 {:?}", gibbs_mat[[0,20]]);
    //floats.extend_from_slice(unsafe{ std::mem::transmute::<&[u8], &[f64]>(&data[0..(nt*8)]) });
    //println!("{:?}", size);
}

#[allow(dead_code)]
pub fn get_ambig(filename: &std::path::Path) -> Vec<u32> {
    let file = File::open(filename).unwrap();
    let buf_reader = BufReader::new(file);

    let mut only_ambig = Vec::<u32>::new();

    for (i, l) in buf_reader.lines().skip(1).enumerate() {
        let s = l.unwrap();
        let mut iter = s.split_ascii_whitespace();
        let u: u32 = iter.next().unwrap().parse().unwrap();
        let a: u32 = iter.next().unwrap().parse().unwrap();
        if u == 0 && a != 0 {
            only_ambig.push(i as u32);
        }
    }
    only_ambig
}

#[allow(dead_code)]
pub fn group_reader(filename: &std::path::Path) -> Vec<Vec<usize>> {
    let file = File::open(filename).unwrap();
    let buf_reader = BufReader::new(file);

    let mut groups = Vec::new();
    for (_i, l) in buf_reader.lines().enumerate() {
        let s = l.unwrap();
        let v: Vec<_> = s.trim().rsplit(',').collect();
        let group: Vec<usize> = v.iter().map(|n| n.parse::<usize>().unwrap()).collect();

        groups.push(group.clone());
    }
    groups
}

// normal variance
/*
fn var_1d(a: ArrayView1<'_, f64>) -> f64 {
    a.var_axis(Axis(0), 0.).into_scalar()
}
*/

// find min and max divide by mean
fn infrv_1d(a: ArrayView1<'_, f64>) -> f64 {
    let mu = a.mean().unwrap() - 1.;
    let var = a.var_axis(Axis(0), 1.).into_scalar();
    //(var) / (mu + 0.1) + 0.01
    if (var - mu) >= 0. {
        (var - mu) / (mu + 5.) + 0.01
    } else {
        0.01
    }
}

// find min and max divide by mean
fn spread1d(a: ArrayView1<'_, f64>) -> f64 {
    let n = a.len() as f64;
    if n == 0. {
        return 0.;
    }
    let mean = a.sum() / n;
    let minimum = a.fold(f64::INFINITY, |m, &x| m.min(x));
    let maximum = a.fold(-f64::INFINITY, |m: f64, &x| m.max(x));
    (maximum - minimum) / mean
}

fn infrv(a: &Array2<f64>, axis: Axis) -> Array1<f64> {
    a.map_axis(axis, infrv_1d)
}

fn spread(a: &Array2<f64>, axis: Axis) -> Array1<f64> {
    a.map_axis(axis, spread1d)
}

pub fn get_infrv_percentile(gibbs_mat: &Array2<f64>, p: f64) -> f64 {
    assert!(0. < p);
    assert!(p < 1.);
    let gibbs_mat_sum = gibbs_mat.sum_axis(Axis(1));
    let gibbs_nz: Vec<_> = gibbs_mat_sum
        .indexed_iter()
        .filter_map(|(index, &item)| if item > 1.0 { Some(index) } else { None })
        .collect();

    let infrv_array = infrv(&gibbs_mat, Axis(1));
    let mut infrv_sort: Vec<f64> = gibbs_nz.iter().map(|i| infrv_array[*i as usize]).collect();
    let n = infrv_sort.len();
    rgsl::sort::vectors::sort(&mut infrv_sort, 1, n);
    rgsl::statistics::quantile_from_sorted_data(&infrv_sort, 1, infrv_sort.len(), p)
}

pub fn endpoints_overdispersed(infrv_array: &Array1<f64>, m: f64, x: usize, y: usize) -> bool {
    infrv_array[x] >= m && infrv_array[y] >= m
}

/*
fn variance(a: &Array2<f64>, axis: Axis) -> Array1<f64> {
    a.map_axis(axis, var_1d)
}
*/

pub fn get_threhold(
    gibbs_mat: &Array2<f64>,
    infrv_quant: f64,
    seed: u64,
    file_list: &FileList,
) -> f64 {
    println!("Calculating threhold");
    let gibbs_mat_sum = gibbs_mat.sum_axis(Axis(1));
    let gibbs_mat_mean = gibbs_mat.mean_axis(Axis(1)).unwrap();
    let gibbs_nz: Vec<_> = gibbs_mat_sum
        .indexed_iter()
        .filter_map(|(index, &item)| if item > 1.0 { Some(index) } else { None })
        .collect();

    let infrv_array = infrv(&gibbs_mat, Axis(1));

    let dat = gibbs_nz
        .iter()
        .map(|i| {
            format!(
                "{}\t{}",
                gibbs_mat_mean[*i as usize], infrv_array[*i as usize]
            )
        })
        .join("\n");
    let infrv_log = file_list.prefix.as_path().join("infrv.log");
    std::fs::write(infrv_log, dat).expect("could not write to the infrv.log");

    let die_roll_log = file_list.prefix.as_path().join("die_roll.log");
    let mut dfile = File::create(die_roll_log).expect("could not create die roll.log");
    // let infrv_array = variance(&gibbs_mat, Axis(1));
    let mut converged = false;
    let starting_num_samples = (gibbs_nz.len() as f64) * 1.;
    println!("\n\nstarting samp : {}\n\n", starting_num_samples);

    let mut starting_num_samples = starting_num_samples as usize;

    let mut old_threshold = 0.0 as f64;
    let mut new_threshold = 0.0 as f64;

    // let mut rng = thread_rng();
    let mut rng = Pcg64::seed_from_u64(seed);
    while !converged {
        //starting_num_samples < gibbs_nz.len(){
        let die_range = Uniform::new(0, gibbs_nz.len());
        let mut roll_die = die_range.sample_iter(&mut rng);

        let mut sampled_infrv = vec![OrderedFloat(0.0); starting_num_samples as usize];
        let mut dice_iter = 0 as usize;
        let mut mean_sum = 0.0f64;

        while dice_iter < starting_num_samples as usize {
            let i1 = roll_die.next().unwrap();
            let mut i2 = roll_die.next().unwrap();
            while i1 == i2 {
                i2 = roll_die.next().unwrap();
            }

            let (t1, t2) = (gibbs_nz[i1], gibbs_nz[i2]);
            if !endpoints_overdispersed(&infrv_array, infrv_quant, t1, t2) {
                continue;
            }
            let s = get_collapse_score(&gibbs_mat, &infrv_array, t1, t2);
            // let s = get_collapse_score(&gibbs_mat, &infrv_array, t1, t2);
            // let s = get_variance_fold_change(&gibbs_mat, &infrv_array, t1, t2);
            // let s = get_infrv_fold_change(&gibbs_mat, &infrv_array, t1, t2);
            let msg = format!(
                "{}\t{}\t{}\t{}\t{}\n",
                gibbs_mat_mean[t1], gibbs_mat_mean[t2], infrv_array[t1], infrv_array[t2], s
            );
            dfile
                .write_all(&msg.into_bytes())
                .expect("could not write to log");
            sampled_infrv[dice_iter] = OrderedFloat(s);
            mean_sum += s;
            dice_iter += 1;

            if dice_iter % 1000 == 1 {
                print!("dice roll: {}\r", dice_iter);
            }
        }
        // calculate threhold
        sampled_infrv.sort();
        let mean = mean_sum / (dice_iter as f64);
        let shifted_samples: Vec<f64> = sampled_infrv
            .iter()
            .map(|s| s.to_f64().unwrap() - mean)
            .collect();
        let shifted_samples_pos: Vec<f64> = shifted_samples
            .into_iter()
            .filter(|&x| x >= 0.0)
            .collect::<Vec<_>>();
        let mid = shifted_samples_pos.len() / 2;
        let median = shifted_samples_pos[mid];
        //let median = sampled_infrv[sampled_infrv.len()/2].to_f64().unwrap();
        new_threshold = mean - (median * 1.48 * 1.95);
        //let sinfrv : Vec<f64> = sampled_infrv.iter().map(|x| x.into_inner()).collect();
        //new_threshold = rgsl::statistics::quantile_from_sorted_data(&sinfrv, 1, sinfrv.len(), 0.025);

        if ((new_threshold - old_threshold) / new_threshold) < 0.001 {
            //- new_threshold).abs() < 1e-3{
            converged = true;
        }

        old_threshold = new_threshold;
        starting_num_samples *= 2;
    }
    new_threshold
}

pub fn get_collapse_score(
    gibbs_mat: &Array2<f64>,
    infrv_array: &Array1<f64>,
    x: usize,
    y: usize,
) -> f64 {
    let infa = infrv_array[x];
    let infb = infrv_array[y];
    let sum = &gibbs_mat.row(x) + &gibbs_mat.row(y);
    //let submat = stack![Axis(0), gibbs_mat.slice(s![x..x+1,..]), gibbs_mat.slice(s![y..y+1,..])];
    //let covmat = submat.cov(1.).unwrap();
    let infsum = infrv_1d(sum.view());
    //let (maxrv, minrv) = if infa < infb { (infb, infa) } else { (infa, infb) };
    infsum - ((infa + infb) * 0.5)
    //infsum - (infa + infb) - covmat[[0,1]]
}

/*
pub fn get_infrv_fold_change(
    gibbs_mat: &Array2<f64>,
    infrv_array: &Array1<f64>,
    x: usize,
    y: usize
) -> f64 {
   let infa = infrv_array[x];
   let infb = infrv_array[y];
   let sum = &gibbs_mat.row(x) + &gibbs_mat.row(y);
   let infsum = infrv_1d(sum.view());
   infsum / (infa + infb)
}
*/

/*
pub fn get_variance_fold_change(
    gibbs_mat: &Array2<f64>,
    variance_array: &Array1<f64>,
    x: usize,
    y: usize
) -> f64 {
   let infa = variance_array[x];
   let infb = variance_array[y];
   let sum = &gibbs_mat.row(x) + &gibbs_mat.row(y);
   let varsum = var_1d(sum.view());
   varsum / (infa + infb + 1.)
}
*/

#[allow(dead_code, clippy::too_many_arguments, clippy::cognitive_complexity)]
pub fn eq_experiment_to_graph(
    exp: &EqClassExperiment,
    gibbs_mat: &mut Array2<f64>,
    eq_class_count: &[u32],
    tolerance: f64,
    thr: f64,
    infrv_quant: f64,
    min_spread: f64,
    delta_file: &mut File,
    unionfind_struct: &mut UnionFind<usize>,
) -> pg::Graph<usize, EdgeInfo, petgraph::Undirected> {
    let start = Instant::now();

    // create a hash of eqclasses
    /*
    println!("Creating hash table for eqclasses");
    let mut trinfo_list = vec![TranscriptInfo::new(); exp.targets.len()];
    let mut golden_set = vec![Vec::<u32>::new(); exp.targets.len()];
    let mut trinfo_map : HashMap<TranscriptInfo, usize> = HashMap::new();
    for (i,x) in exp.classes.iter().enumerate() {
        if i % 10000 == 1 {
            print!("{}\r",i);
            io::stdout().flush().unwrap();
        }
        let ns = x.0;
        let ws = x.1;

        let wsrounded : Vec<i32> = ws.iter().map(|&w| (w * 1000f32).round() as i32).collect();
        for j in 0..ns.len(){
            let a = ns[j];
            let w = wsrounded[j];
            trinfo_list[a as usize].eqlist.push(i);
            trinfo_list[a as usize].weights.push(w);
        }
    }
    // insert into the hash
    let mut golden_set_index = 0 as usize ;
    for (i,tinfo) in trinfo_list.iter().enumerate(){
        if i == 107880 || i == 241729{
            println!("{:?}", tinfo);
        }
        if tinfo.eqlist.len() == 0 { continue; }
        match trinfo_map.get(&tinfo) {
            Some(&stored_tinfo_index) => {
                golden_set[stored_tinfo_index].push(i as u32);
            },
            None => {
                trinfo_map.insert(tinfo.clone(),golden_set_index);
                golden_set[golden_set_index].push(i as u32);
                golden_set_index += 1;
            }
        }
    }

    for (_,tlist_bk) in golden_set.iter().enumerate(){
        let mut tlist = tlist_bk.clone();
        tlist.sort();
        if tlist.len() > 1{
            let source = tlist[0];
            for t in tlist.iter().skip(1){
                let to_add = gibbs_mat.index_axis(Axis(0), *t as usize).to_owned();
                let mut s = gibbs_mat.slice_mut(s![source as usize,..]);
                s += &to_add;
                unionfind_struct.union(source as usize, *t as usize);
            }
        }
    }
    */

    // code for refinery
    println!("Creating partition refinery");
    let mut part_cache = HashSet::new();
    let part_start = Instant::now();
    let mut valid_transcripts = vec![false; exp.targets.len()];
    let mut part = Partition::simple(exp.targets.len());
    for (i, x) in exp.classes.iter().enumerate() {
        if i % 10000 == 1 {
            print!("{}\r", i);
            io::stdout().flush().unwrap();
        }
        let ns = x.0;
        let ws = x.1;

        if ns.len() < 2 {
            part.refine(&[ns[0] as usize]);
            continue;
        }

        //let thresh = 0.1 * (1.0 / ns.len() as f32);
        //let retained : std::vec::Vec<usize> = (0..ns.len()).filter_map(
        //    |j| if ws[j as usize] >= thresh { Some(j as usize) } else { None} ).collect();

        let _wsrounded: Vec<i32> = ws.iter().map(|&w| (w * 1000f32).round() as i32).collect();
        let mut pair_vec = Vec::with_capacity(ns.len());
        for j in 0..ns.len() {
            pair_vec.push((ns[j], OrderedFloat(ws[j])));
            valid_transcripts[ns[j] as usize] = true;
        }
        pair_vec.sort_by_key(|k| k.1);

        //let mut partition_sets = Vec::new();

        let mut j = 0 as usize;
        let mut tmp_vec = Vec::with_capacity(ns.len());

        while j < pair_vec.len() - 1 {
            let mut diff = pair_vec[j + 1].1.to_f64().unwrap() - pair_vec[j].1.to_f64().unwrap();
            tmp_vec.clear();
            tmp_vec.push(pair_vec[j].0 as usize);
            while diff < tolerance {
                tmp_vec.push(pair_vec[j + 1].0 as usize);
                valid_transcripts[pair_vec[j].0 as usize] = true;
                valid_transcripts[pair_vec[j + 1].0 as usize] = true;
                j += 1;
                if j < pair_vec.len() - 1 {
                    diff = pair_vec[j + 1].1.to_f64().unwrap() - pair_vec[j].1.to_f64().unwrap();
                } else {
                    break;
                }
            }
            //partition_sets.push(tmp_vec.clone());
            if part_cache.insert(tmp_vec.clone()) {
                part.refine(&tmp_vec[..]);
            }
            j += 1;
        }
        if j < pair_vec.len() && part_cache.insert(vec![pair_vec[j].0 as usize]) {
            part.refine(&[pair_vec[j].0 as usize]);
            //valid_transcripts[pair_vec[j].0 as usize] = true;
            //partition_sets.push(vec![pair_vec[j].0 as usize]);
        }

        //let mut tmp_vec = Vec::with_capacity(ns.len());
        //let mut last_idx = 0 as usize;
        //let mut ce = pair_vec[last_idx];
        //loop {
        //   tmp_vec.clear();
        //   tmp_vec.extend(pair_vec.iter().skip(last_idx).take_while(|x| (x.1.to_f64().unwrap() - ce.1.to_f64().unwrap()) < 1e-3).map(|x| x.0 as usize));
        //   part.refine(&tmp_vec[..]);
        //   last_idx += tmp_vec.len();
        //   if last_idx >= pair_vec.len() { break; }
        //   ce = pair_vec[last_idx];
        //}
        //println!("{:?}",partition_sets);
    }

    let part_vec = part.iter().collect::<Vec<_>>();
    let mut golden_collapses = 0;
    for (_, p) in part_vec.iter().enumerate() {
        if p.len() > 1 {
            //println!("{:?}", p);
            if valid_transcripts[p[0]] {
                if p.len() > 10 {
                    println!("{},{}", p.len(), p[0]);
                }
                let mut tlist = p.to_vec();
                tlist.sort();
                let source = tlist[0];
                for t in tlist.iter().skip(1) {
                    let to_add = gibbs_mat.index_axis(Axis(0), *t as usize).to_owned();
                    let mut s = gibbs_mat.slice_mut(s![source as usize, ..]);
                    s += &to_add;
                    unionfind_struct.union(source as usize, *t as usize);
                    golden_collapses += 1;
                }
            } else if p.len() > 10 {
                println!("{}", p.len());
            }
        }
    }
    println!("Number of golden collapses {}", golden_collapses);
    println!("The refinery code ran for {:?}", part_start.elapsed());

    let mut og = pg::Graph::<usize, EdgeInfo, petgraph::Undirected>::new_undirected();
    for (i, _n) in exp.targets.iter().enumerate() {
        let idx = og.add_node(i);
        // the index assigned by the graph should be the
        // order in which we add these
        debug_assert_eq!(i, idx.index());
    }

    // compute the mean of the elemenrs
    let gibbs_mat_mean = gibbs_mat.mean_axis(Axis(1)).unwrap();
    let gibbs_mat_spread = spread(&gibbs_mat, Axis(1));

    // apply threashold
    let filtered_indices_spread: Vec<u32> = gibbs_mat_spread
        .indexed_iter()
        .filter_map(|(index, &item)| {
            if item > min_spread && gibbs_mat_mean[index] > 1.0 {
                Some(index as u32)
            } else {
                None
            }
        })
        .collect();

    let filtered_indices_exp: Vec<u32> = gibbs_mat_spread
        .indexed_iter()
        .filter_map(|(index, &_item)| {
            if gibbs_mat_mean[index] > 1.0 {
                Some(index as u32)
            } else {
                None
            }
        })
        .collect();

    let mut filtered_indices_vec = vec![false; exp.targets.len()];
    let mut filtered_indices_mean_vec = vec![false; exp.targets.len()];

    for i in filtered_indices_spread {
        filtered_indices_vec[i as usize] = true;
    }
    for i in filtered_indices_exp {
        filtered_indices_mean_vec[i as usize] = true;
    }

    //let shape = gibbs_mat.shape() ;
    let infrv_array = infrv(&gibbs_mat, Axis(1));
    //let infrv_array = variance(&gibbs_mat, Axis(1));

    //let start_corr = Instant::now();
    // let gibbs_filtered_mat = keep_rows(&gibbs_mat, &filtered_indices);
    // println!("Duration for computing sub matrix {:?}", start_corr.elapsed());
    //let gibbs_corr_mat = gibbs_filtered_mat.pearson_correlation().unwrap();
    // println!("Computing correlation {:?}", start_corr.elapsed());
    //correlation(data1: &[f64], stride1: usize, data2: &[f64], stride2: usize, n: usize)

    for (i, x) in exp.classes.iter().enumerate() {
        let ns = x.0;
        let ws = x.1;
        let eq_count = x.2;
        assert!(eq_count == eq_class_count[i]);

        let thresh = 0.1 * (1.0 / ns.len() as f32);

        let retained: std::vec::Vec<usize> = (0..ns.len())
            .filter_map(|j| {
                if ws[j as usize] >= thresh {
                    Some(ns[j] as usize)
                } else {
                    None
                }
            })
            .collect();

        for a in 0..retained.len() {
            let mut na = retained[a] as usize;
            let na_root = unionfind_struct.find(na);
            if na_root != na {
                na = na_root;
            }

            for nb in retained.iter().skip(a + 1) {
                let mut nbd = *nb as usize;
                let nb_root = unionfind_struct.find(nbd);
                if nb_root != nbd {
                    nbd = nb_root;
                }

                if na == nbd {
                    continue;
                }

                if endpoints_overdispersed(&infrv_array, infrv_quant, na, nbd)
                    && (filtered_indices_vec[na] || filtered_indices_vec[nbd])
                //&&
                //  (filtered_indices_mean_vec[na] && filtered_indices_mean_vec[*nb])
                {
                    let (u, v) = if na < nbd { (na, nbd) } else { (nbd, na) };

                    let va = pg::graph::NodeIndex::new(u);
                    let vb = pg::graph::NodeIndex::new(v);
                    let e = og.find_edge(va, vb);
                    match e {
                        Some(ei) => {
                            let mut ew = og.edge_weight_mut(ei).unwrap();
                            ew.count += eq_count;
                            ew.eqlist.push(i);
                        }
                        None => {
                            // only add the edge if the correlation is sufficientl
                            // small
                            let delta = get_collapse_score(&gibbs_mat, &infrv_array, na, *nb);
                            // let delta = get_variance_fold_change(&gibbs_mat, &infrv_array, na, *nb);
                            //let delta = get_infrv_fold_change(&gibbs_mat, &infrv_array, na, *nb);

                            let s = format!("{}\t{}\t{}\n", na, *nb, delta);
                            delta_file
                                .write_all(s.as_bytes())
                                .expect("failed to write to delta.log file");
                            //    if (u == 116212) && (v == 116212){
                            //        println!("=================={}================", delta);
                            //    }
                            if delta < thr {
                                og.add_edge(
                                    va,
                                    vb,
                                    EdgeInfo {
                                        infrv_gain: delta,
                                        count: eq_count,
                                        state: -1,
                                        eqlist: vec![i],
                                    },
                                );
                            }
                        }
                    }
                } // if not filtered and count big enough
            }
        }
        //println!("{:?} :: {:?} :: [{:?}]", x.0, x.1, x.2);
    }

    let og2 = og.filter_map(
        |_ni, n| Some(*n),
        |ei, e| {
            let (a, b) = og.edge_endpoints(ei).unwrap();
            let min_mean = gibbs_mat_mean[a.index()].min(gibbs_mat_mean[b.index()]);
            if (e.count as f64) > min_mean {
                Some(EdgeInfo {
                    infrv_gain: e.infrv_gain,
                    count: e.count,
                    state: -1,
                    eqlist: e.eqlist.clone(),
                })
            } else {
                None
            }
        },
    );
    println!("Prev node count: {}", og.node_count());
    println!("Prev edge count: {}", og.edge_count());
    println!("New node count: {}", og2.node_count());
    println!("New edge count: {}", og2.edge_count());
    //println!("# cc : {}", pg::algo::connected_components(&og));
    println!("Elapsed time for computing graph {:?}", start.elapsed());
    og2
}

// find intersection
fn intersect<T: std::cmp::PartialOrd + std::cmp::Ord + Copy>(a: &[T], b: &[T]) -> std::vec::Vec<T> {
    let mut ia = a.iter();
    let mut ib = b.iter();
    let mut va = ia.next();
    let mut vb = ib.next();
    let mut out = std::vec::Vec::with_capacity(std::cmp::min(a.len(), b.len()));
    while va.is_some() && vb.is_some() {
        if *(va.unwrap()) < *(vb.unwrap()) {
            va = ia.next();
        } else {
            let not_less_or_equal = match (*vb.unwrap()).partial_cmp(&*(va.unwrap())) {
                None | Some(Ordering::Equal) | Some(Ordering::Greater) => true,
                _ => false,
            };
            if not_less_or_equal {
                out.push(*(va.unwrap()));
                va = ia.next();
            }
            vb = ib.next();
        }
    }
    out.dedup();
    out
}

#[allow(dead_code)]
fn union<T: std::cmp::PartialOrd + Copy>(a: &[T], b: &[T]) -> std::vec::Vec<T> {
    let mut ia = a.iter();
    let mut ib = b.iter();
    let mut va = ia.next();
    let mut vb = ib.next();
    let mut out = std::vec::Vec::with_capacity(a.len() + b.len());
    while va.is_some() || vb.is_some() {
        if va.is_some() && vb.is_some() {
            if *(va.unwrap()) < *(vb.unwrap()) {
                out.push(*(va.unwrap()));
                va = ia.next();
            } else {
                let not_less_or_equal = match (*vb.unwrap()).partial_cmp(&*(va.unwrap())) {
                    None | Some(Ordering::Greater) => true,
                    _ => false,
                };
                if not_less_or_equal {
                    out.push(*(va.unwrap()));
                    va = ia.next();
                    vb = ib.next();
                } else {
                    out.push(*(vb.unwrap()));
                    vb = ib.next();
                }
            }
        } else if va.is_some() {
            out.push(*(va.unwrap()));
            va = ia.next();
        } else {
            out.push(*(vb.unwrap()));
            vb = ib.next();
        }
    }
    out
}

#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord)]
struct EdgeState {
    infrv_gain: OrderedFloat<f64>,
    source: usize,
    target: usize,
    state: i32,
}

pub fn verify_graph(
    eq_class_count: &[u32],
    og: &mut pg::Graph<usize, EdgeInfo, petgraph::Undirected>,
) {
    for e in og.edge_references() {
        let eweight = e.weight();
        let c = eweight.count;
        let mut vc = 0usize;
        for eq in &eweight.eqlist {
            vc += eq_class_count[*eq] as usize;
        }
        if vc != c as usize {
            println!("verified edge count = {}, but stored count = {}", vc, c);
        }
    }
}

//util::work_on_component(&eq_class, &gibbs_array, &og, &comp);
#[allow(clippy::too_many_arguments, clippy::cognitive_complexity)]
pub fn work_on_component(
    eq_class_count: &[u32],
    gibbs_mat: &mut Array2<f64>,
    gibbs_mat_mean: &mut Array1<f64>,
    unionfind_struct: &mut UnionFind<usize>,
    og: &mut pg::Graph<usize, EdgeInfo, petgraph::Undirected>,
    comp: &[pg::graph::NodeIndex],
    num_collapses: &mut usize,
    thr: f64,
    infrv_quant: f64,
    cfile: &mut File,
) {
    // make a set of edges to be visited
    let mut infrv_array = infrv(&gibbs_mat, Axis(1));
    //let mut infrv_array = variance(&gibbs_mat, Axis(1));
    //let shape = gibbs_mat.shape().to_vec() ;

    let mut heap = BinaryHeap::new_min();

    for node in comp.iter() {
        for e in og.edges(*node) {
            // let (mut min_node, mut max_node) = (node_deref, next_node);
            let source = e.source().index();
            let target = e.target().index();
            let w = e.weight().infrv_gain;
            let (u, v) = if source > target {
                (target, source)
            } else {
                (source, target)
            };
            assert!(u < v, "source = {}, target = {}", u, v);

            heap.push(EdgeState {
                infrv_gain: OrderedFloat(w),
                source: u,
                target: v,
                state: -1,
            });
        }
    }

    while let Some(EdgeState {
        infrv_gain,
        source,
        target,
        state,
    }) = heap.pop()
    {
        // take the next available edge that can be collapsed
        // correlation less than threashold
        if infrv_gain < OrderedFloat(thr) {
            // println!("mincorr {}\r",corr);
            // if the edge count satisfies the criteria
            let source_node = pg::graph::NodeIndex::new(source);
            let target_node = pg::graph::NodeIndex::new(target);

            assert!(
                source_node.index() < target_node.index(),
                "source = {}, source index = {}, target = {}, target index = {}",
                source,
                source_node.index(),
                target,
                target_node.index()
            );

            let u_to_v_id_opt = og.find_edge(source_node, target_node);

            // only proceed if this edge still exists
            if u_to_v_id_opt.is_none() {
                continue;
            }

            let u_to_v_id = u_to_v_id_opt.unwrap();
            let u_to_v_info = og.edge_weight_mut(u_to_v_id).unwrap();

            // only proceed if the state of this edge when
            // placed on the heap is still the current state
            if u_to_v_info.state != state {
                continue;
            }

            // the min posterior mean among u and v
            let min_mean = gibbs_mat_mean[source].min(gibbs_mat_mean[target]);

            // the count along this edge must be large enough to "matter"
            if f64::from(u_to_v_info.count) >= min_mean {
                let msg = format!(
                    "{}\t{}\t{}\t{}\t{}\n",
                    source, target, infrv_array[source], infrv_array[target], infrv_gain
                );
                //println!("{}",msg);
                cfile
                    .write_all(&msg.into_bytes())
                    .expect("could not write into collapse log");

                // one collapse guaranteed
                *num_collapses += 1;

                // remove the collapsed edge first so that
                // u is not in neighbors(v) and
                // v is not in neighbors(u)
                og.remove_edge(u_to_v_id);

                // it's candidate for collapse
                // when we collapse we modify the
                // following
                // i. gibbs_mat
                // ii. gibbs_mat_mean
                // iii. Each neighbor of u and v
                // iv. heap
                // v. unionfind_array
                unionfind_struct.union(source, target);

                let to_add = gibbs_mat.index_axis(Axis(0), target).to_owned();
                let mut s = gibbs_mat.slice_mut(s![source, ..]);
                s += &to_add;
                infrv_array[source] = infrv_1d(s.view());
                gibbs_mat_mean[source] = s.sum() / (s.len() as f64);

                // update correlation for (u*v) to new and existing neighbors
                let mut source_adj: Vec<usize> = og
                    .neighbors(source_node)
                    .map(|n| n.clone().index())
                    .collect();
                let mut target_adj: Vec<usize> = og
                    .neighbors(target_node)
                    .map(|n| n.clone().index())
                    .collect();

                source_adj.sort();
                target_adj.sort();

                source_adj.dedup();
                target_adj.dedup();
                let common_ids = intersect(&source_adj, &target_adj);

                // edge u-v is collapsed
                // if u-x is already an edge, but v-x is *not*
                // then we don't update count
                for x in source_adj.iter() {
                    if common_ids.contains(x) {
                        continue;
                    }
                    // it is only neighbor of u
                    // so we need to update correlation

                    let xn = pg::graph::NodeIndex::new(*x);
                    let u_to_x_inner = og.find_edge(source_node, xn).unwrap();
                    let mut u_to_x_info_inner = og.edge_weight_mut(u_to_x_inner).unwrap();
                    let curr_state = u_to_x_info_inner.state;

                    let delta = get_collapse_score(&gibbs_mat, &infrv_array, source, *x);
                    // let delta = get_variance_fold_change(&gibbs_mat, &infrv_array, source, *x);
                    // let delta = get_infrv_fold_change(&gibbs_mat, &infrv_array, source, *x);

                    u_to_x_info_inner.infrv_gain = delta;
                    u_to_x_info_inner.state += 1;

                    // update heap
                    if delta < thr && endpoints_overdispersed(&infrv_array, infrv_quant, source, *x)
                    {
                        let (a, b) = if source > *x {
                            (*x, source)
                        } else {
                            (source, *x)
                        };
                        assert!(
                            a != b,
                            "1. source = {}, target = {}, source = {}, target = {}",
                            a,
                            b,
                            source,
                            target
                        );
                        heap.push(EdgeState {
                            infrv_gain: OrderedFloat(delta),
                            source: a,
                            target: b,
                            state: curr_state + 1,
                        });
                    }
                }

                // if there is v -- x but not
                // u -- x then we copy the properties
                // of v -- x into our new u*v -- x
                for x in target_adj.iter() {
                    if common_ids.contains(x) {
                        continue;
                    }

                    let xn = pg::graph::NodeIndex::new(*x);
                    let v_to_x_inner = og.find_edge(target_node, xn).unwrap();
                    let v_to_x_info_inner = og.edge_weight(v_to_x_inner).unwrap();

                    let v_to_x_count = v_to_x_info_inner.count;
                    let v_to_x_eqlist = &v_to_x_info_inner.eqlist.to_vec();

                    let delta = get_collapse_score(&gibbs_mat, &infrv_array, source, *x);
                    // let delta = get_variance_fold_change(&gibbs_mat, &infrv_array, source, *x);
                    // let delta = get_infrv_fold_change(&gibbs_mat, &infrv_array, source, *x);

                    let new_state = -1 as i32;

                    og.add_edge(
                        source_node,
                        xn,
                        EdgeInfo {
                            infrv_gain: delta,
                            count: v_to_x_count,
                            state: new_state,
                            eqlist: v_to_x_eqlist.to_vec(),
                        },
                    );

                    // update heap
                    if delta < thr && endpoints_overdispersed(&infrv_array, infrv_quant, source, *x)
                    {
                        let (a, b) = if source > *x {
                            (*x, source)
                        } else {
                            (source, *x)
                        };
                        assert!(
                            a != b,
                            "2. source = {}, target = {}, source = {}, target = {}",
                            a,
                            b,
                            source,
                            target
                        );
                        heap.push(EdgeState {
                            infrv_gain: OrderedFloat(delta),
                            source: a,
                            target: b,
                            state: new_state,
                        });
                    }
                    og.remove_edge(v_to_x_inner);
                }

                // it's a triangle
                for x in common_ids.iter() {
                    let xn = pg::graph::NodeIndex::new(*x);

                    // for the u - x edge
                    let u_to_x_inner = og.find_edge(source_node, xn).unwrap();
                    // for the v - x edge
                    let v_to_x_inner = og.find_edge(target_node, xn).unwrap();

                    let v_to_x_count: u32;
                    let v_to_x_eq: Vec<usize>;
                    {
                        let v_to_x_info = og.edge_weight(v_to_x_inner).unwrap();
                        v_to_x_count = v_to_x_info.count;
                        v_to_x_eq = v_to_x_info.eqlist.clone();
                    }

                    let mut u_to_x_info = og.edge_weight_mut(u_to_x_inner).unwrap();

                    // v_to_x_eq.sort();
                    let intersecting_eqlist = intersect(&v_to_x_eq, &u_to_x_info.eqlist);
                    let curr_state = u_to_x_info.state;

                    let mut sum = 0 as u32;
                    for i in intersecting_eqlist.iter() {
                        sum += eq_class_count[*i as usize];
                    }
                    let tot_current_count = v_to_x_count + u_to_x_info.count;

                    // TODO(@hiraksarkar) : once confident, we can remove this check
                    if sum > tot_current_count {
                        println!("sum = {}, tot_current_count = {}", sum, tot_current_count);
                        println!(
                            "u-x : {:?}, v-x : {:?}, intersection : {:?}",
                            u_to_x_info.eqlist, v_to_x_eq, intersecting_eqlist
                        );
                        std::process::exit(1);
                    }

                    // FIXME(@hiraksarkar) : bootleg intersection --- replace with the proper function
                    // once we have it
                    //for c in &v_to_x_eq {
                    //    if !intersecting_eqlist.contains(c) {
                    //        u_to_x_info.eqlist.push(*c);
                    //    }
                    //}
                    //u_to_x_info.eqlist.sort();
                    u_to_x_info.eqlist = union(&u_to_x_info.eqlist.to_vec(), &v_to_x_eq);

                    let final_count = tot_current_count - sum;

                    let delta = get_collapse_score(&gibbs_mat, &infrv_array, source, *x);
                    // let delta = get_variance_fold_change(&gibbs_mat, &infrv_array, source, *x);
                    // let delta = get_infrv_fold_change(&gibbs_mat, &infrv_array, source, *x);

                    u_to_x_info.infrv_gain = delta;
                    u_to_x_info.count = final_count;
                    u_to_x_info.state += 1;

                    // update heap
                    if delta < thr && endpoints_overdispersed(&infrv_array, infrv_quant, source, *x)
                    {
                        let (a, b) = if source > *x {
                            (*x, source)
                        } else {
                            (source, *x)
                        };
                        assert!(
                            a != b,
                            "3. a = {}, b = {}, source = {}, target = {}",
                            a,
                            b,
                            source,
                            target
                        );
                        heap.push(EdgeState {
                            infrv_gain: OrderedFloat(delta),
                            source: a,
                            target: b,
                            state: curr_state + 1,
                        });
                    }
                    og.remove_edge(v_to_x_inner);
                }
            }
        }
    }
    //println!("curr state {}", curr_state);
}

pub fn parse_eq(filename: &std::path::Path) -> Result<EqClassExperiment, io::Error> {
    //let start = Instant::now();
    let file = File::open(filename).expect("equivalence class file does not exist");
    let reader: Box<dyn Read> = if filename.ends_with("eq_classes.txt.gz") {
        Box::new(GzDecoder::new(file))
    } else {
        Box::new(file)
    };
    let mut buf_reader = BufReader::new(reader);
    let mut buf = String::new();

    let mut exp = EqClassExperiment::new();

    buf_reader
        .read_line(&mut buf)
        .expect("Cannot read first line");
    buf.pop();
    let num_target: usize = buf.parse().unwrap();
    exp.ntarget = num_target;
    buf.clear();

    buf_reader
        .read_line(&mut buf)
        .expect("Cannot read second line");
    buf.pop();
    let num_eq: usize = buf.parse().unwrap();
    //let count: u64 = buf.parse().unwrap();

    exp.neq = num_eq;

    let mut tnames = Vec::<String>::with_capacity(num_target);

    for _ in 0..num_target {
        buf.clear();
        buf_reader
            .read_line(&mut buf)
            .expect("could read target name");
        buf.pop();
        tnames.push(buf.to_string());
    }

    exp.targets = tnames;

    //let mut pb = pbr::ProgressBar::new(count);
    //pb.format("╢▌▌░╟");

    for _ in 0..num_eq {
        buf.clear();
        buf_reader
            .read_line(&mut buf)
            .expect("could read eq. class");
        buf.pop();
        let mut iter = buf.split_ascii_whitespace();
        let nt: usize = iter.next().unwrap().parse().unwrap();
        let mut tv = Vec::<u32>::with_capacity(nt);
        let mut wv = Vec::<f32>::with_capacity(nt);
        for _ in 0..nt {
            tv.push(iter.next().unwrap().parse().unwrap());
        }
        for _ in 0..nt {
            wv.push(iter.next().unwrap().parse().unwrap());
        }
        let c: u32 = iter.next().unwrap().parse().unwrap();
        exp.add_class(&mut tv, &mut wv, c);
        //pb.inc();
    }
    //pb.finish_print("done");
    //let duration = start.elapsed();
    Ok(exp)
    // make graph from these
}


pub fn matrix_reader(
    input: &str,
    num_cells: usize,
    num_genes: usize,
    expressions: &mut Vec<Vec<u8>>,
    bit_vecs: &mut Vec<Vec<u8>>,
    tier_fraction_vec: &mut Vec<f32>,
) -> Result<bool, io::Error> {
    println!("Using {} as input EDS file\n", input);
    println!(
        "Using {} Rows (cells) and {} Columns (features)",
        num_cells, num_genes
    );

    let num_bit_vecs: usize = round::ceil(num_genes as f64 / 8.0, 0) as usize;
    let mut total_molecules = 0;
    let mut total_exp_values = 0;

    {
        let mut count = 0;
        let file_handle = File::open(input)?;
        let mut file = GzDecoder::new(file_handle);

        for _ in 0..num_cells {
            let mut bit_vec = vec![0; num_bit_vecs];
            file.read_exact(&mut bit_vec[..])?;
            let mut num_ones = 0;
            for bits in bit_vec.iter() {
                num_ones += bits.count_ones();
            }
            bit_vecs.push(bit_vec);

            //let mut expression: Vec<u8> = vec![0; 4 * (num_ones as usize)];
            //let mut float_buffer: Vec<f32> = vec![0.0_f32; num_ones as usize];
            //file.read_exact(&mut expression[..])?;

            let mut expression: Vec<u8> = vec![0; 1 * (num_ones as usize)];

            file.read_exact(&mut expression[..])?;

            let cell_count: u32 = expression.iter().map(|&x| x as u32).sum();
            total_molecules += cell_count;

            //expr.slice_mut(s![i as usize, ..])
            //    .assign(&Array::from(expression));
            expressions.push(expression);

            count += 1;
            total_exp_values += num_ones;
            if count % 100 == 0 {
                print!("\r Done Reading {} cells", count);
                io::stdout().flush()?;
            }
        }
    }

    println!("\n");
    assert!(
        expressions.len() == num_cells,
        "rows and quants file size mismatch"
    );

    println!("Found Total {:.2} molecules", total_molecules);
    println!("Found Total {:.2} expressed entries", total_exp_values);
    println!(
        "w/ {:.2} Molecules/cell",
        total_molecules as f32 / num_cells as f32
    );
    
    for (cell_id, exp) in expressions.into_iter().enumerate() {
        let bit_vec = &bit_vecs[cell_id];
        let mut fids: Vec<usize> = Vec::new();

        for (feature_id, flag) in bit_vec.into_iter().enumerate() {
            if *flag != 0 {
                for (offset, j) in format!("{:8b}", flag).chars().enumerate() {
                    match j {
                        '1' => fids.push((8 * feature_id) + offset),
                        _ => (),
                    };
                }
            }
        }

        assert!(
            fids.len() == exp.len(),
            format!("#positions {} doesn't match with #expressed features {}",
                    fids.len(), exp.len())
        );
        // mtx_data = format!("cell{}", cell_id + 1);
        let mut zero_counter = 0;
        for (index, count) in exp.into_iter().enumerate() {
            assert!(
                fids[index] < num_genes,
                format!("{} position > {}", fids[index], num_genes)
            );

            while zero_counter != fids[index] {
                zero_counter += 1;
                // mtx_data.push_str(&format!(",0"));
            }

            zero_counter += 1;
            if *count == 3u8 {
                tier_fraction_vec[fids[index]] += 1.0;
            }
            // mtx_data.push_str(&format!(",{}", count));
        }

        while zero_counter < num_genes {
            zero_counter += 1;
            // mtx_data.push_str(&format!(",0"));
        }

        // mtx_data.push_str(&format!("\n"));
        // file.write_all(mtx_data.as_bytes())?;
    }

    // for i in &mut tier_fraction_vec {
    //     *i /= num_cells as f32;
    // }
    tier_fraction_vec.iter_mut().for_each(|i| *i /= num_cells as f32);
    Ok(true)
}

pub fn parse_bfh(
    alevin_info: &AlevinMetaData, 
    t2g_file_name: &std::path::Path,
    _tier_mat: & Vec<Vec<u8>>,
) -> Result<BFHEqClassExperiment, io::Error> {

    // Get transcript to gene mapping
    let mut t2gmap = HashMap::new();
    let t2g_file = File::open(t2g_file_name).expect("transcript to gene mapping file does not exist");
    #[derive(Deserialize)]
    struct T2GPair {
        transcript_name: String,
        gene_name: String,
    };

    let mut rdr = csv::ReaderBuilder::new()
        .has_headers(false)
        .delimiter(b'\t')
        .from_reader(t2g_file);
    for result in rdr.deserialize() {
        // Notice that we need to provide a type hint for automatic
        // deserialization.
        let record: T2GPair = result?;
        t2gmap.insert(record.transcript_name,record.gene_name);
    }

    // --- reading trnascripts, cell barcodes and number of 
    // --- equivalence classes
    // read bfh file
    let bfh_file = File::open(alevin_info.bfh_file.clone()).expect("can not open bfh file");
    let mut buf_reader = BufReader::new(bfh_file);
    let mut buf = String::new();
    let mut exp = BFHEqClassExperiment::new();
    // num_transcripts
    buf_reader
        .read_line(&mut buf)
        .expect("Cannot read first line");
    buf.pop();
    let num_targets : usize = buf.parse().unwrap();
    buf.clear();
    // num_cells
    buf_reader
        .read_line(&mut buf)
        .expect("Cannot read second line");
    buf.pop();
    let num_cells : usize = buf.parse().unwrap();
    buf.clear();
    // num_equivalence classes
    buf_reader
        .read_line(&mut buf)
        .expect("Cannot read second line");
    buf.pop();
    let neq : usize = buf.parse().unwrap();
    buf.clear();

    // convert the transcripts to gene maps
    let mut tnames = Vec::<String>::with_capacity(num_targets);
    for _ in 0..num_targets {
        buf.clear();
        buf_reader
            .read_line(&mut buf)
            .expect("could read target name");
        buf.pop();
        tnames.push(buf.to_string());
    }
    // reads cell names
    let mut cnames = Vec::<String>::with_capacity(num_cells);
    for _ in 0..num_cells {
        buf.clear();
        buf_reader
            .read_line(&mut buf)
            .expect("could read target name");
        buf.pop();
        cnames.push(buf.to_string());
    }

    exp.neq = neq;

    for j in 0..neq {
        buf.clear();
        buf_reader
            .read_line(&mut buf)
            .expect("could read eq. class");
        buf.pop();
        let mut iter = buf.split_ascii_whitespace();
        let num_labels : usize = iter.next().unwrap().parse().unwrap();
        // things to push to class
        let mut gene_labels = Vec::<usize>::new() ;
        let mut cell_ids = Vec::<usize>::new() ;
        // let mut tiers = Vec::<u16>::new() ;       
        // read the tids
        for _ in 0..num_labels {
            let tid : usize = iter.next().unwrap().parse().unwrap();
            // transcript name
            let tname = tnames[tid].clone() ;
            // gene name
            if let Some(gname) = t2gmap.get(&tname) {
                if let Some(gid) = alevin_info.feature_map.get(gname) {
                    // insert this gene id
                    gene_labels.push(*gid);
                }
            }
            // may contain duplicates
            gene_labels.sort_unstable();
            gene_labels.dedup();
        }

        // tot_num_reads
        let tot_num_reads : u32 = iter.next().unwrap().parse().unwrap();
        // num_bc
        let num_bcs : usize = iter.next().unwrap().parse().unwrap();
        for _ in 0..num_bcs {
            // cell barcode
            let bc : usize = iter.next().unwrap().parse().unwrap();
            let bc_name = cnames[bc].clone();
            if let Some(cell_id) = alevin_info.cell_barcode_map.get(&bc_name) {
                cell_ids.push(*cell_id);
            }
            // number of umi,count pair
            let num_umi : usize = iter.next().unwrap().parse().unwrap();
            for _ in 0..num_umi {
                let _umi_seq = iter.next();
                let _umi_cnt = iter.next();
            }
        }
        // add the class
        exp.add_class(&mut gene_labels,&mut cell_ids, tot_num_reads);
        if j % 100 == 0 {
            print!("\r Done Reading {} equivalence classes", j);
            io::stdout().flush()?;
        }
    }

    Ok(exp)
}