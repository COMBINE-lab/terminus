use serde::{Deserialize, Serialize};
use std::fs::*;
use std::hash::Hash;
use std::io::BufReader;
use std::path::PathBuf;

#[derive(Debug)]
pub struct ConsensusFileList {
    pub cons_nwk_file: PathBuf,
    pub cluster_bp_splits_file: PathBuf,
    pub merged_groups_file: PathBuf,
    //pub groups_length: PathBuf,
}
impl ConsensusFileList {
    pub fn new(dname: String) -> ConsensusFileList {
        let dir = PathBuf::from(dname);
        if !dir.as_path().exists() {
            panic!("The directory {} did not exist", dir.to_str().unwrap());
        }
        if !dir.as_path().is_dir() {
            panic!(
                "The path {} did not point to a valid directory",
                dir.to_str().unwrap()
            );
        }
        
        ConsensusFileList {
            cluster_bp_splits_file: dir.as_path().join("cluster_bipart_splits.txt"),
            cons_nwk_file: dir.as_path().join("cluster_nwk.txt"),
            merged_groups_file: dir.as_path().join("merged_groups_length.txt"),
      //      groups_length: dir.as_path().join("groups_length.txt")
        }
    }
}

#[derive(Debug)]
pub struct FileList {
    pub prefix: PathBuf,
    pub quant_file: PathBuf,
    pub bootstrap_file: PathBuf,
    pub ambig_file: PathBuf,
    pub mi_file: PathBuf,
    pub eq_file: PathBuf,
    pub names_tsv_file: PathBuf,
    pub cmd_file: PathBuf,
    pub collapsed_log_file: PathBuf,
    pub golden_collapses_log_file: PathBuf,
    pub allele_collapses_log_file: PathBuf,
    pub group_file: PathBuf,
    pub collapse_order_file: PathBuf,
    pub delta_file: PathBuf,
    pub cluster_file: PathBuf,
    pub gene_cluster_file: PathBuf,
    pub group_bp_splits_file: PathBuf,
    pub cluster_bp_splits_file: PathBuf,
    pub group_nwk_file: PathBuf,
    pub mgroup_nwk_file: PathBuf,
    pub param_log_file: PathBuf,
    
}

// construct the files
impl FileList {
    pub fn new(dname: String) -> FileList {
        let dir = PathBuf::from(dname);
        if !dir.as_path().exists() {
            panic!("The directory {} did not exist", dir.to_str().unwrap());
        }
        if !dir.as_path().is_dir() {
            panic!(
                "The path {} did not point to a valid directory",
                dir.to_str().unwrap()
            );
        }
        let parent = dir.as_path().join("aux_info");

        let mut eq_name = "eq_classes.txt";
        let mi_path = parent.join("meta_info.json");
        if mi_path.exists() {
            let file = File::open(mi_path);
            let reader = BufReader::new(file.unwrap());
            let jd: MetaInfo = serde_json::from_reader(reader).unwrap();

            eq_name = if jd.eq_class_properties.contains(&"gzipped".to_string()) {
                "eq_classes.txt.gz"
            } else {
                "eq_classes.txt"
            };
        }

        FileList {
            prefix: dir.clone(),
            ambig_file: parent.join("ambig_info.tsv"),
            mi_file: parent.join("meta_info.json"),
            quant_file: dir.as_path().join("quant.sf"),
            eq_file: parent.join(eq_name),
            bootstrap_file: parent.join("bootstrap").join("bootstraps.gz"),
            names_tsv_file: parent.join("bootstrap").join("names.tsv.gz"),
            cmd_file: dir.as_path().join("cmd_info.json"),
            cluster_file: dir.as_path().join("clusters.txt"),
            collapsed_log_file: dir.as_path().join("collapsed.log"),
            golden_collapses_log_file: dir.as_path().join("golden_collapses.log"),
            allele_collapses_log_file: dir.as_path().join("allele_collapses.log"),
            group_file: dir.as_path().join("groups.txt"),
            //group_order_file: dir.as_path().join("order.txt"),
            collapse_order_file: dir.as_path().join("collapse_order.json"),
            delta_file: dir.as_path().join("delta.log"),
            gene_cluster_file: dir.as_path().join("gene_cluster.log"),
            group_bp_splits_file: dir.as_path().join("group_bipart_splits.txt"),
            cluster_bp_splits_file: dir.as_path().join("cluster_bipart_splits.txt"),
            group_nwk_file: dir.as_path().join("group_nwk.txt"),
            mgroup_nwk_file: dir.as_path().join("mgroup_nwk.txt"),
            param_log_file: dir.as_path().join("param_log_file.json"),
            //cluster_nwk_file: dir.as_path().join("cluster_nwk.txt"),
        }
    }
}

#[derive(PartialEq, Eq, Hash, Clone, Debug, Default)]
pub struct TranscriptInfo {
    pub eqlist: Vec<usize>,
    pub weights: Vec<i32>,
}

impl TranscriptInfo {
    pub fn new() -> TranscriptInfo {
        TranscriptInfo {
            eqlist: Vec::<usize>::new(),
            weights: Vec::<i32>::new(),
        }
    }
}

pub struct EdgeInfo {
    pub infrv_gain: f64,
    pub count: u32,
    pub state: i32,
    pub eqlist: Vec<usize>,
}

#[derive(Debug, Default, Clone)]
pub struct EqList {
    pub offsets: Vec<usize>,
    pub labels: Vec<u32>,
    pub weights: Vec<f32>,
    pub counts: Vec<u32>,
}

pub struct IterEqList<'a> {
    inner: &'a EqList,
    pos: usize,
}

/// Iterator over an EqList that yields
/// each (label, weight, count) tuple in turn.
///
impl<'a> Iterator for IterEqList<'a> {
    type Item = (&'a [u32], &'a [f32], u32);
    fn next(&mut self) -> Option<Self::Item> {
        let i = self.pos;
        if (i + 1) >= self.inner.offsets.len() {
            return None;
        }
        self.pos += 1;
        let p = self.inner.offsets[i];
        let l = self.inner.offsets[(i + 1)] - p;
        Some((
            &self.inner.labels[p..(p + l)],
            &self.inner.weights[p..(p + l)],
            self.inner.counts[i],
        ))
    }
}

impl EqList {
    pub fn iter(&self) -> IterEqList {
        IterEqList {
            inner: self,
            pos: 0,
        }
    }

    pub fn add_class(&mut self, labs: &mut Vec<u32>, weights: &mut Vec<f32>, count: u32) {
        let len = weights.len();
        self.offsets.push(self.offsets.last().unwrap() + len);
        self.labels.append(labs);
        self.weights.append(weights);
        self.counts.push(count);
    }

    pub fn new() -> EqList {
        EqList {
            offsets: vec![0_usize],
            labels: Vec::<u32>::new(),
            weights: Vec::<f32>::new(),
            counts: Vec::<u32>::new(),
        }
    }
}

#[derive(Debug, Clone)]
pub struct EqClassExperiment {
    pub targets: Vec<String>,
    pub ntarget: usize,
    pub neq: usize,
    pub classes: EqList,
}

impl EqClassExperiment {
    /// Add an equivalence class to the set of equivalence classes for this experiment
    ///
    /// # Arguments
    ///
    /// * `labs` - A vector of the transcript ids for this eq class
    /// * `weights` - A vector of the conditional probability weights for the
    ///               transcripts in `labs`.
    /// * `count` - The number of fragments that are equivalent with respect to this class.
    ///
    pub fn add_class(&mut self, labs: &mut Vec<u32>, weights: &mut Vec<f32>, count: u32) {
        self.classes.add_class(labs, weights, count);
    }

    pub fn new() -> EqClassExperiment {
        EqClassExperiment {
            targets: Vec::<String>::new(),
            ntarget: 0,
            neq: 0,
            classes: EqList::new(),
        }
    }
}

impl Default for EqClassExperiment {
    fn default() -> Self {
        Self::new()
    }
}

#[derive(Debug, Deserialize, Clone)]
#[allow(non_snake_case)]
pub struct TxpRecord {
    pub Name: String,
    pub Length: u32,
    pub EffectiveLength: f32,
    pub TPM: f32,
    pub NumReads: f32,
}

// Iterator over an EqList that yields
/// each (label, weight, count) tuple in turn.
///
//impl<'a> Iterator for IterTxpRecord<'a> {
//    type Item = (&'a String, &'a u32, &'a f32, &'a f32, &'a f32);
//    fn next(&mut self) -> Option<Self::Item> {
//        let i = self.pos;
//        if i >= self.inner.0.len() {
//            return None;
//        }
//        self.pos += 1;
//        Some((
//            &self.inner.Name,
//            self.inner.Length,
//            self.inner.EffectiveLength,
//            self.inner.TPM,
//            self.inner.NumReads,
//        ))
//    }
//}

impl TxpRecord {
    // pub fn iter(&self) -> IterTxpRecord {
    //     IterTxpRecord {
    //         inner: self,
    //         pos: 0,
    //     }
    // }

    pub fn assign(&mut self, other: &TxpRecord) {
        self.Name = other.Name.clone();
        self.Length = other.Length;
        self.EffectiveLength = other.EffectiveLength;
        self.TPM = other.TPM;
        self.NumReads = other.NumReads;
    }

    pub fn new() -> TxpRecord {
        TxpRecord {
            Name: "".to_string(),
            Length: 0u32,
            EffectiveLength: 0.0f32,
            TPM: 0.0f32,
            NumReads: 0.0f32,
        }
    }
}

impl Default for TxpRecord {
    fn default() -> Self {
        Self::new()
    }
}

#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct TerminusInfo {
    pub num_nontrivial_groups: usize,
    pub num_targets_in_nontrivial_group: usize,
}

#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct MetaInfo {
    pub num_valid_targets: u32,
    pub serialized_eq_classes: bool,
    pub num_bootstraps: u32,
    pub num_eq_classes: u32,
    pub eq_class_properties: Vec<String>,
    pub samp_type: String,
}

// pub struct MetaInfo {
//     pub salmon_version : String,
//     pub samp_type: String,
//     pub opt_type: String,
//     pub quant_errors: Vec<String>,
//     pub num_libraries: u32,
//     pub library_types: Vec<String>,
//     pub frag_dist_length: u32,
//     pub seq_bias_correct: bool,
//     pub gc_bias_correct: bool,
//     pub num_bias_bins: u32,
//     pub mapping_type: String,
//     pub num_valid_targets : u32,
//     pub num_decoy_targets : u32,
//     pub num_eq_classes: u32,
//     pub serialized_eq_classes: bool,
//     pub eq_class_properties: Vec<String>,
//     pub length_classes: Vec<u32>,
//     pub index_seq_hash: String,
//     pub index_name_hash: String,
//     pub index_seq_hash512: String,
//     pub index_name_hash512: String,
//     pub index_decoy_seq_hash: String,
//     pub index_decoy_name_hash: String,
//     pub num_bootstraps: u32,
//     pub num_processed: u32,
//     pub num_mapped: u32,
//     pub num_decoy_fragments: u32,
//     pub num_dovetail_fragments: u32,
//     pub num_fragments_filtered_vm: u32,
//     pub num_alignments_below_threshold_for_mapped_fragments_vm: usize,
//     pub percent_mapped: f32,                                                                                                                                                                                                                                                                                          "call": "quant",
//     pub start_time: String,
//     pub end_time: String,
// }

pub struct TxpRecordSet {
    pub samples: Vec<Vec<TxpRecord>>,
}
