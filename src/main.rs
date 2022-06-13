pub mod salmon_types;
mod util;
mod collapse;
pub mod binary_tree;

extern crate serde_stacker;
extern crate serde_json;
extern crate serde_pickle;


use std::collections::{HashMap, HashSet, BTreeMap};
use std::fs::*;
use std::io::Write;
use std::io::{self, BufRead, BufReader};
use std::process;

use clap::{App, AppSettings, Arg, ArgMatches, SubCommand};
use ndarray::prelude::*;
use num_format::{Locale, ToFormattedString};
use petgraph as pg;
use petgraph::algo::{connected_components, tarjan_scc};
use petgraph::unionfind::UnionFind;
use rayon::prelude::*;
use std::path::PathBuf;


use serde::Deserialize;
use serde_json::Value;


// Name of the program, to be used in diagnostic messages.
static PROGRAM_NAME: &str = "terminus";

/// Exit the program, printing an error message on stderr, and returning
/// a specific error code. The program name is prefixed onto the front of
/// the error message. If logging is enabled we also write the error
/// message to the log file.
#[allow(unused)]
fn exit_with_error(status: i32, message: &str) {
    //error!(target: "log_messages", "{}", message);
    writeln!(
        &mut std::io::stderr(),
        "{} ERROR: {}!",
        PROGRAM_NAME,
        message
    )
    .unwrap();
    std::process::exit(status);
}

fn do_group(sub_m: &ArgMatches) -> Result<bool, io::Error> {
    //let mut groups: Vec<Vec<usize>> = Vec::new();
    //let dir_paths : Vec<_> = sub_m.values_of("dirs").unwrap().collect();

    let dname: String = sub_m.value_of("dir").unwrap().to_string();
    let prefix: String = sub_m.value_of("out").unwrap().to_string();
    create_dir_all(prefix.clone())?;

    // let mut unionfind_vec = Vec::<UnionFind<_>>::with_capacity(dir_paths.len()) ;
    // let mut gibbs_array_vec = Vec::<Array2<f64>>::with_capacity(dir_paths.len());
    // let mut meta_info_array = Vec::<salmon_types::MetaInfo>::with_capacity(dir_paths.len());
    // let mut num_global_targrts = 0u32 ;

    let seed = sub_m
        .value_of("seed")
        .unwrap()
        .parse::<u64>()
        .expect("generate random values from the seed");
    let min_spread = sub_m
        .value_of("min-spread")
        .unwrap()
        .parse::<f64>()
        .expect("could not convert min-spread to float value");

    let tolerance = sub_m
        .value_of("tolerance")
        .unwrap()
        .parse::<f64>()
        .expect("could not convert tolerance to float value");
    
    let thr_bool = sub_m
        .value_of("thresh")
        .unwrap()
        .parse::<bool>()
        .expect("could not parse thr_bool");

    println!("------input configuration------");
    println!("seed : {}", seed);
    println!("min-spread : {}", min_spread);
    println!("tolerance : {}", tolerance);
    println!("dir : {}", dname);
    let compo: Vec<&str> = dname.rsplit('/').collect();
    //println!("{:?}",compo);
    let experiment_name = compo[0];
    let mut prefix_path = prefix;
    prefix_path.push('/');
    prefix_path.push_str(experiment_name);
    let file_list = salmon_types::FileList::new(dname.to_string());
    
    // create output directory
    
    println!("output folder: {}", prefix_path);
    println!("------------------------------");
    // create
    create_dir_all(prefix_path.clone())?;
    let file_list_out = salmon_types::FileList::new(prefix_path);

    // Load the gibbs samples
    let x = util::parse_json(&file_list.mi_file).unwrap();

    let mut gibbs_array =
        Array2::<f64>::zeros((x.num_valid_targets as usize, x.num_bootstraps as usize));
    util::read_gibbs_array(&file_list.bootstrap_file, &x, &mut gibbs_array);
    let mut gibbs_mat_mean = gibbs_array.mean_axis(Axis(1)).unwrap();

    // if a2g exists also dumps gene level groups
    let mut allele2txp = PathBuf::from(sub_m.value_of("a2t").unwrap().to_string());
    let asemode:bool = allele2txp.as_path().is_file();
    if asemode {
        println!("Alleles would be collapsed according to the file: {:?}",allele2txp.to_str().unwrap());
    }

    // if t2g exists restrict equivalence classes to gene level groups
    let mut transcript2gene = PathBuf::from(sub_m.value_of("t2g").unwrap().to_string());
    let txpmode:bool = transcript2gene.as_path().is_file();
    if txpmode {
        println!("Txps within a gene would be collapsed using : {:?}", transcript2gene.to_str().unwrap());
    }

    // take the transcript to gene mapping
    // this will also create a map from transcript id
    // to gene id

    println!("parsing eqfile {:?}", file_list.eq_file);
    let eq_class = util::parse_eq(&file_list.eq_file).unwrap();
    println!("length of eqclass {:?}", eq_class.neq);
    let mut eq_class_counts = vec![0_u32; eq_class.neq];
    // let mut i = 0_usize;
    for (i, eq) in eq_class.classes.iter().enumerate() {
        eq_class_counts[i] = eq.2;
    }
    let mut collapse_order: Vec<binary_tree::TreeNode> = Vec::new();
    for i in 0..eq_class.ntarget {
        collapse_order.push(binary_tree::TreeNode::create_leaf(i.to_string()));
    }

    // fill targets from eq_class
    let mut tnames:HashMap<String, usize> = HashMap::new();
    for i in 0..eq_class.ntarget {
        tnames.insert(eq_class.targets[i].clone(), i);
    }
    let ntarget = eq_class.ntarget;

    let mut gene2allele_map:HashMap<usize, Vec<usize>> = HashMap::new();
    let mut allele2gene_map:Vec<usize> = Vec::new();

    let mut txp2allele_map:HashMap<usize, Vec<usize>> = HashMap::new();
    let mut allele2txp_map:Vec<usize> = Vec::new();
    
    if asemode {
        util::get_map_bw_ent(&mut allele2txp, &mut txp2allele_map, &mut allele2txp_map, &tnames);
        let nalleles = allele2txp_map.len();
        if nalleles != ntarget {
            panic!("number of alleles {} not equal to number of txps in eq class file {}", nalleles, ntarget);
        }
    }

    if txpmode {
        util::get_map_bw_ent(&mut transcript2gene, &mut gene2allele_map, &mut allele2gene_map, &tnames);    
        let ntxps = allele2gene_map.len();
        if ntxps != ntarget {
            panic!("number of txps {} not equal to number of txps in eq class file {}", ntxps, ntarget);
        }
    }

    if asemode {
        let nalleles = allele2txp_map.len();
        if txpmode {
            let ntxps = allele2gene_map.len();
            if nalleles != ntxps {
                panic!("number of alleles {} not equal to number of txps in eq class file {}", nalleles, ntxps);
            }
        }
    }    
        

    let inf_perc = 0.25f64;
    let p = util::get_infrv_percentile(&gibbs_array, inf_perc);
    //let p = 2.48675518f64;
    println!("the {}% of infRV was : {}", inf_perc * 100., p);
    let mut thr = 1e7 as f64;
    if thr_bool {
        thr = util::get_threhold(&gibbs_array, p, seed, &file_list_out);
    }
    
   // let thr = -12.958284980475035f64;
    // thr = thr * 0.75;
    // thr = 0.645;
    println!("threshold: {}", thr);
    //let dpath = Path::new(file_list_out.delta_file.clone());
    let mut dfile =
        File::create(file_list_out.delta_file.clone()).expect("could not create collapse.log");
    let mut unionfind_struct = UnionFind::new(eq_class.ntarget);
    let mut group_order:Vec<String>=Vec::with_capacity(eq_class.ntarget);
    for i in 0..eq_class.ntarget {
        group_order.push(i.to_string())
    }
    
   
    let mut ed_vec = util::eq_experiment_to_graph(
        &eq_class,
        &mut gibbs_array,
        &eq_class_counts,
        tolerance,
        thr,
        p,
        min_spread,
        &mut dfile,
        &mut unionfind_struct,
        &allele2gene_map,
        &txp2allele_map,
        &mut group_order,
        &mut collapse_order
    );
    


    // connected coponents
//     let num_connected_components = connected_components(&gr);
//     println!("#Connected components {:?}", num_connected_components);

    let mut num_collapses = 0_usize;

//     //let cpath = Path::new(file_list_out.collapsed_log_file.clone());
    let mut cfile = File::create(file_list_out.collapsed_log_file.clone())
        .expect("could not create collapse.log");

//     let gcomp: Vec<petgraph::prelude::NodeIndex> = gr
//         .node_indices()
//         .map(|x| petgraph::graph::NodeIndex::new(x.index()))
//         .collect();
    util::work_on_component(
        &eq_class_counts,
        &mut gibbs_array,
        &mut gibbs_mat_mean,
        &mut unionfind_struct,
        &mut ed_vec,
        &mut num_collapses,
        thr,
        p,
        &mut cfile,
        &mut group_order,
        &mut collapse_order
    );

    // //write down the groups
    let mut groups = HashMap::new();
    //let mut grouped_set = HashSet::new();
    for i in 0..(x.num_valid_targets as usize) {
        let root = unionfind_struct.find(i);
        if root != i {
            groups.entry(root).or_insert_with(Vec::new).push(i);
        }
    }

    // println!("Number of collapsed transcripts from conn components with 2 {}", num_collapses_2.to_formatted_string(&Locale::en));
    println!(
        "Number of collapses {}",
        num_collapses.to_formatted_string(&Locale::en)
    );
    //let _res = util::write_modified_quants(&groups, &grouped_set, &file_list_out, &gibbs_array, &x, &rec, &collapsed_dim);
    
    let mut gfile = File::create(file_list_out.group_file).expect("could not create groups.txt");
    let mut co_file = File::create(file_list_out.collapse_order_file).expect("could not create collapse order file");
    let mut nwk_file = File::create(file_list_out.group_nwk_file).expect("could not create group order file");
    // let mut co = File::create("co_file.txt").expect("could not create collapse order file");
    let _write = util::group_writer(&mut gfile, &groups);
   // let _write = util::order_group_writer(&mut gofile, &group_order, &groups);
    let _write = util::collapse_order_writer(&mut co_file, &mut nwk_file, &groups, &collapse_order);
    //let _write = util::co_id_writer(&mut co, &groups, &collapse_order);
    
    // let file = File::open(ff);
    // let reader = BufReader::new(file.unwrap());
    //let dd: HashMap<usize,binary_tree::TreeNode> = serde_pickle::from_reader(reader).unwrap();
        
    Ok(true)
}

fn do_collapse(sub_m: &ArgMatches) -> Result<bool, io::Error> { 

    let num_threads = sub_m
        .value_of("threads")
        .unwrap()
        .parse::<usize>()
        .expect("could not parse --threads option");
    // set the number of rayon threads
    rayon::ThreadPoolBuilder::new()
        .num_threads(num_threads)
        .build_global()
        .unwrap_or_else(|_| panic!("could not set number of rayon threads to {}", num_threads));

    //let dir_paths: Vec<&str> = sub_m.values_of("dirs").unwrap().collect();
    let sal_dir: String = sub_m.value_of("dirs").unwrap().to_string();
    let md = metadata(sal_dir.clone()).unwrap_or_else(|_| panic!("Invalid directory {}", sal_dir));
    
    let mut sal_dir_paths = read_dir(sal_dir)?
        .map(|res| res.map(|e| e.path()))
        .filter(|res| res.as_ref().unwrap().is_dir())
        .collect::<Result<Vec<_>, io::Error>>()?;
    sal_dir_paths.sort();
    println!("{:?}", sal_dir_paths);
    let mut dir_paths: Vec<&str> = Vec::new();
    for entry in sal_dir_paths.iter() {
        dir_paths.push(entry.as_path().to_str().unwrap());
    }
    
    
    let prefix: String = sub_m.value_of("out").unwrap().to_string();

    
    //let mut bipart_counter: HashMap<String, u32> = HashMap::new();
    let mut bipart_counter: HashMap<String, HashMap<String, u32>> = HashMap::new();
    let mut global_graph = pg::Graph::<usize, u32, petgraph::Undirected>::new_undirected();
    let mut group_keys:Vec<String> = Vec::new();
    let mut l = 0; // num of groups in 1st 
    let mut ntxps = 0; // num of transcripts
    let mut tnames:Vec<String> = Vec::new();
    // add edges
    for (i, dname) in dir_paths.iter().enumerate() {
        //let dname = dname.as_path().to_str().unwrap();
        let compo: Vec<&str> = dname.rsplit('/').collect();
        if compo[0].chars().next().unwrap() == '.' {
            continue;
        }
        let experiment_name = compo[0];
        println!("experiment name {}", experiment_name);
        
        //let mut bipart_counter: HashMap<String, u32> = HashMap::new();
        let mut dir_bipart_counter: HashMap<String, HashMap<String, u32>> = HashMap::new(); // Storing counts of each bipartition
        //let mut group_bipart: HashMap<String, Vec<String>> = HashMap::new(); // Storing all bipartitions per group
        
        let mut prefix_path = prefix.clone();
        prefix_path.push('/');
        prefix_path.push_str(experiment_name);
        create_dir_all(prefix_path.clone())?;
        let file_list_out = salmon_types::FileList::new(prefix_path);
        let file = File::open(file_list_out.collapse_order_file);
        let reader = BufReader::new(file.unwrap());
    
        let mut deserializer = serde_json::Deserializer::from_reader(reader);
        deserializer.disable_recursion_limit();

    //println!("{:?}", deserializer);
        let collapse_order: HashMap<String,binary_tree::TreeNode> = HashMap::deserialize(&mut deserializer).unwrap(); // can be replaced by vector
        for (key, node) in &collapse_order  {
            let node = collapse_order.get(key).unwrap();
            let req_group = binary_tree::sort_group_id(&node.id);
            //let node_vec = group_bipart.entry(node.id.clone()).or_insert(Vec::<String>::new());
            let dir_group_key = dir_bipart_counter.entry(req_group.clone()).or_insert(HashMap::new());
            let overall_group_key = bipart_counter.entry(req_group.clone()).or_insert(HashMap::new());
            let node_set:HashSet<u32> = node.id.clone()
                                        .split("_")
                                        .map(|x| x.parse::<u32>()
                                        .unwrap()).collect();
            //binary_tree::compute_bipart_count(node, &mut bipart_counter, &mut dir_bipart_counter, &node_set, node_vec);
            group_keys.push(req_group.clone());
            binary_tree::compute_bipart_count2(node, overall_group_key, dir_group_key);
        }
        if i == 0 {
            l = dir_bipart_counter.len();
            let file_old = salmon_types::FileList::new(dname.to_string());
            let eq_class = util::parse_eq(&file_old.eq_file).unwrap();
            ntxps = eq_class.ntarget;
         
            // for i in 0..eq_class.targets.len()
            // {
            //     tnames.push(i.to_string());
            // }
            tnames.extend(eq_class.targets.clone());
        }
        println!("Number of groups in {} are {}", dname, dir_bipart_counter.len());
        let mut bipart_file = File::create(file_list_out.group_bp_splits_file).expect("could not create group bp splits");
        let _f = util::mapTrait::bipart_writer(&dir_bipart_counter, &mut bipart_file, &tnames);

    } // all files
    println!("Total number of groups are {}", bipart_counter.len());
    let m_groups = sub_m
    .value_of("merge_groups")
    .unwrap()
    .parse::<bool>()
    .expect("could not parse --merge_groups option");

    let m_type = sub_m
    .value_of("merge_type")
    .unwrap()
    .to_string();
    

    if m_groups {
        println!("{}",m_type);
        if m_type == "BP"{
            bipart_counter = collapse::merge_groups(&bipart_counter, ntxps);
            println!("Groups after merging {:?}", bipart_counter.len());        
            let mut co = File::create("co_file.txt").expect("could not create collapse order file");
            let _write = util::mapTrait::bipart_writer(&bipart_counter, &mut co, &tnames);
        }
        else{
            if vec![String::from("phylip"), String::from("stelar")].contains(&m_type) {
                let all_groups:Vec<String> = bipart_counter.keys().cloned().collect();
                collapse::use_phylip(&dir_paths, &prefix, &all_groups, ntxps, &m_type);
            }
        }
    }
    
    // filter based on the threshold
    let consensus_thresh = sub_m
        .value_of("consensus-thresh")
        .unwrap()
        .parse::<f64>()
        .expect("could not parse --consensus-thresh option");
    if m_type == "BP" {
        let half_length = (dir_paths.len() as f64 * consensus_thresh).floor() as u32;
        let global_filtered_graph = global_graph.filter_map(
            |_, n| Some(*n),
            |_, &e| {
                if (e as usize) >= half_length as usize{
                    Some(e)
                } else {
                    None
                }
            },
        );
        for (_,group) in bipart_counter.iter_mut(){
            //let group_bpart = bipart_counter.get_mut(&group.clone()).unwrap();
            group.retain(|key, value| {
                *value >= half_length 
            });
        }
        let mut total_in_group = 0;
        let mut num_group = 0;
        for (bipart_split, _) in bipart_counter.iter() {
            total_in_group += bipart_split.split("_").collect::<Vec<&str>>().len() + 1; //+1 for bp
            num_group += 1;
        }
        let term_info = salmon_types::TerminusInfo {
            num_nontrivial_groups: num_group,
            num_targets_in_nontrivial_group: total_in_group,
        };
        println!(
            "total number of transcripts in a non-trivial group : {}",
            total_in_group
        );
        println!("total consensus clusters : {}", num_group);
        let file_list_out = salmon_types::ConsensusFileList::new(prefix.clone());
        let mut bipart_file = File::create(file_list_out.cluster_bp_splits_file.clone()).expect("could not create cluster bp splits");
        let _f = util::mapTrait::bipart_writer(&bipart_counter, &mut bipart_file, &tnames);
    }

        
    // if t2g exists also dumps gene level groups
    // let transcript2gene = PathBuf::from(sub_m.value_of("t2g").unwrap().to_string());
    // println!(
    //     "transcript2gene file: {:?}",
    //     transcript2gene.to_str().unwrap()
    // );
    // if transcript2gene.as_path().is_file() {
    //     // get name of the transcripts from any directory
    //     println!("=============Reducing to gene groups=============");
    //     let mut t2gmap: HashMap<String, String> = HashMap::new();
    //     let mut genemap: HashMap<String, u32> = HashMap::new();
    //     let genenames = util::get_t2g(&transcript2gene, &mut genemap, &mut t2gmap);
    //     println!("{} genes exist in the file", genenames.len());

    //     let file_list = salmon_types::FileList::new((dir_paths[0]).to_string());
    //     let x = util::parse_json(&file_list.mi_file).expect("json file could not be parsed");
    //     let rec =
    //         util::parse_quant(&file_list.quant_file, &x).expect("quant file could not be parsed");

    //     let mut genevec: Vec<u32> = vec![0u32; rec.len()];
    //     let mut genevecpresent: Vec<bool> = vec![false; rec.len()];
    //     let mut notfound = 0;
    //     for i in 0..rec.len() {
    //         let tname = rec[i].Name.clone();
    //         // println!("Searching for {:?}", tname);
    //         match t2gmap.get(&tname) {
    //             Some(gname) => match genemap.get(gname) {
    //                 Some(geneid) => {
    //                     genevec[i] = *geneid;
    //                     genevecpresent[i] = true;
    //                 }
    //                 None => {
    //                     println!("Not found {:?}, {:?}", tname, gname);
    //                 }
    //             },
    //             None => {
    //                 //println!("transcript name not found {}", tname);
    //                 notfound += 1;
    //             }
    //         }
    //     }
    //     if notfound > 0 {
    //         println!("{} transcripts not in {:?}", notfound, transcript2gene);
    //     }
    //     let mut global_gene_graph = pg::Graph::<usize, u32, petgraph::Undirected>::new_undirected();
    //     if global_gene_graph.node_count() == 0 {
    //         for i in 0..genenames.len() {
    //             let idx = global_gene_graph.add_node(i as usize);
    //             // the index assigned by the graph should be the
    //             // order in which we add these
    //             debug_assert_eq!(i as usize, idx.index());
    //         }
    //     }
    //     println!("Initial graph constructed");
    //     for edge in global_filtered_graph.raw_edges() {
    //         let source = edge.source().index();
    //         let end = edge.target().index();
    //         if !(genevecpresent[source] && genevecpresent[end]) {
    //             continue;
    //         }
    //         let na = genevec[source];
    //         let nb = genevec[end];
    //         let va = pg::graph::NodeIndex::new(na as usize);
    //         let vb = pg::graph::NodeIndex::new(nb as usize);
    //         let e = global_gene_graph.find_edge(va, vb);
    //         match e {
    //             Some(ei) => {
    //                 let ew = global_gene_graph
    //                     .edge_weight_mut(ei)
    //                     .expect("edge weight not found");
    //                 *ew += 1;
    //             }
    //             None => {
    //                 global_gene_graph.add_edge(va, vb, 1);
    //             }
    //         }
    //     }
    //     // components
    //     let mut comps_gene: Vec<Vec<_>> = tarjan_scc(&global_gene_graph);
    //     // comps_gene.sort_by(|v, w| v.len().cmp(&w.len()));
    //     comps_gene.sort_by_key(|v| v.len());
    //     println!("Done reducing to gene level groups");
    //     // let mut gfile = File::create(file_list.gene_cluster_file).expect("could not create groups.txt");
    //     // let _write = util::gene_writer(&mut gfile, &comps_gene, &genenames);
    //     dir_paths.clone().into_par_iter().for_each(|dname| {
    //         let compo: Vec<&str> = dname.rsplit('/').collect();
    //         let experiment_name = compo[0];
    //         let mut prefix_path = prefix.clone();
    //         prefix_path.push('/');
    //         prefix_path.push_str(experiment_name);
    //         let file_list_out = salmon_types::FileList::new(prefix_path);
    //         let mut gfile =
    //             File::create(file_list_out.gene_cluster_file).expect("could not create groups.txt");
    //         let _write = util::gene_writer(&mut gfile, &comps_gene, &genenames);
    //     });
    // }

    // // components
    // let mut comps: Vec<Vec<_>> = tarjan_scc(&global_filtered_graph);
    // comps.sort_by(|v, w| v.len().cmp(&w.len()));
    // comps.sort_by_key(|v| v.len());
    // // write a json file containing
    // // new number of nodes, number of components, etc
    
    //let mut l = term_info;
    // // for each experiment call the writer
    // println!("=============Writing collapsed output=============");
    // dir_paths.into_par_iter().for_each(|dname| {
    //     // getting prepared for output
    //     let compo: Vec<&str> = dname.rsplit('/').collect();
    //     let experiment_name = compo[0];
    //     let mut prefix_path = prefix.clone();
    //     prefix_path.push('/');
    //     prefix_path.push_str(experiment_name);
    //     let file_list_out = salmon_types::FileList::new(prefix_path.to_string());
    //     let file_list = salmon_types::FileList::new(dname.to_string());
    //     // create output directory
    //     println!("Sample {}", experiment_name);

    //     let bpath = std::path::Path::new(&prefix_path);
    //     let term_info_path = bpath.join("terminus_info.json");
    //     let estr = format!("could not write terminus_info.json in {}", prefix_path);
    //     ::serde_json::to_writer(&File::create(term_info_path).expect(&estr), &term_info)
    //         .expect(&estr);
    //     // load old files
    //     // let x = util::parse_json(&file_list.mi_file).unwrap();
    //     // copy cmd_info.json from old location to new
    //     std::fs::copy(&file_list.cmd_file, &file_list_out.cmd_file)
    //         .expect("Could not copy cmd_info.json.");
    //     let x = util::parse_json(&file_list.mi_file).unwrap();
    //     let rec = util::parse_quant(&file_list.quant_file, &x).unwrap();
    //     let mut gibbs_array =
    //         Array2::<f64>::zeros((x.num_valid_targets as usize, x.num_bootstraps as usize));
    //     util::read_gibbs_array(&file_list.bootstrap_file, &x, &mut gibbs_array);
        
    //     let mut bipart_file = File::create(file_list_out.cluster_bp_splits_file.clone()).expect("could not create cluster bp splits");
    //     let _f = util::mapTrait::bipart_writer(&bipart_counter, &mut bipart_file, &tnames);
    // // //     //call the writer
    //     let _res =
    //         util::write_quants_from_components(&comps, &file_list_out, &gibbs_array, &x, &rec);
    // });
    
    Ok(true)
}

// The entry point of the program
// that reads from the files and build
// graphs or later produces the collapsed
// files.

fn main() -> io::Result<()> {
    let matches = App::new("Terminus")
	.setting(AppSettings::ArgRequiredElseHelp)
        .version("0.1.53")
        .author("Sarkar et al.")
        .about("Data-driven grouping of transcripts to reduce inferential uncertainty")
        .subcommand(
            SubCommand::with_name("group")
            .about("perform per-sample grouping of transcripts; required prior to consensus collapse.")
            .arg(
                Arg::with_name("dir")
                    .long("dir")
                    .short("d")
                    .required(true)
                    .takes_value(true)
                    .help("directory to read input from")
            )
            .arg(
                Arg::with_name("min-spread")
                    .long("min-spread")
                    .short("m")
                    .takes_value(true)
                    .default_value("0.1")
                    .help("the minimum spread a transcript must exhibit to enable an \
                          attached edge to be a collapse candidate")
            )
            .arg(
                Arg::with_name("tolerance")
                    .long("tolerance")
                    .takes_value(true)
                    .default_value("0.001")
                    .help("The allowable difference between the weights of transcripts \
                          in same equivalence classes to treat them as identical")
            )
            .arg(
                Arg::with_name("seed")
                    .long("seed")
                    .takes_value(true)
                    .default_value("10")
                    .help("The allowable difference between the weights of transcripts \
                          in same equivalence classes to treat them as identical")
            )
            .arg(
                Arg::with_name("a2t")
                    .long("a2t")
                    .takes_value(true)
                    .default_value("")
                    .help("Mapping allele to transcript")
            )
            .arg(
                Arg::with_name("t2g")
                    .long("t2g")
                    .takes_value(true)
                    .default_value("")
                    .help("Mapping transcript to gene")
            )
            .arg(
                Arg::with_name("thresh")
                    .long("thr")
                    .takes_value(true)
                    .default_value("true")
                    .help("do we use a threshold for collapsing")
            )
            .arg(
                Arg::with_name("out")
                    .long("out")
                    .short("o")
                    .required(true)
                    .takes_value(true)
                    .requires("dir")
                    .help("prefix where output would be written")
            )
        )
        .subcommand(
            SubCommand::with_name("collapse")
            .about("analyze a collection of per-sample groups, and produce a consensus grouping.")
            .arg(
                Arg::with_name("dirs")
                    .long("dirs")
                    .short("d")
                    .required(true)
               //     .multiple(true)
                    .takes_value(true)
                    .help("direcotories to read the group files from")
            )
            .arg(
                Arg::with_name("out")
                    .long("out")
                    .short("o")
                    .required(true)
                    .takes_value(true)
                    .requires("dirs")
                    .help("prefix where output would be written")
            )
            .arg(
                Arg::with_name("t2g")
                    .long("t2g")
                    .takes_value(true)
                    .requires("dirs")
                    .default_value("")
                    .help("use this to group to genes")
            )
            .arg(
                Arg::with_name("threads")
                    .long("threads")
                    .short("t")
                    .takes_value(true)
                    .default_value("8")
                    .help("number of threads to use when writing down the collapsed quantifications")
            )
            .arg(
                Arg::with_name("consensus-thresh")
                .long("consensus-thresh")
                .short("c")
                .takes_value(true)
                .default_value("0.5")
                .help("threshold for edge consensus")
            )
            .arg(
                Arg::with_name("merge_groups")
                .long("merge-groups")
                .short("m")
                .takes_value(true)
                .default_value("true")
                .help("Merge groups")
            )
            .arg(
                Arg::with_name("merge_type")
                .long("merge_type")
                .takes_value(true)
                .default_value("phylip")
                .help("Merging method")
            )
        ).get_matches();

    pretty_env_logger::init_timed();

    match matches.subcommand() {
        ("group", Some(sub_m)) => {
            do_group(&sub_m).expect("Grouping failed");
        }
        ("collapse", Some(sub_m)) => {
            do_collapse(&sub_m).expect("Grouping failed");
        }
        _ => unreachable!(),
    }

    Ok(())
}
