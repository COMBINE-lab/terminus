use serde::{Serialize, Deserialize};
use petgraph::unionfind::UnionFind;
use std::collections::{HashMap, HashSet, BTreeMap};
use std::iter::FromIterator;
use crate::binary_tree::{sort_group_id};

pub fn create_union_find(g:&[String], ntxps:usize) -> UnionFind<usize> {
    let mut unionfind_struct = UnionFind::new(ntxps);
    let mut visited:Vec<i32> = vec![-1; ntxps];
    let mut count = 0;
    for (i,group) in g.iter().enumerate() {
        let mut g_set:Vec<usize> = group.clone().split("_").map(|x| x.parse::<usize>().unwrap()).collect();
        let source = g_set[0];
        if visited[source as usize] == -1 {
            count += 1;
            visited[source as usize] = 0
        }
        for t in g_set.iter().skip(1){
            unionfind_struct.union(source, *t);
            if visited[*t as usize] == -1 {
                count += 1;
                visited[*t as usize] = 0
            }
        }
    }
    if count != ntxps {
        panic!("The number of expected transcripts {} do not match the counted transcripts from the groups{}", ntxps, count);
    }
    return unionfind_struct;
}

pub fn find_group_inds(groups:&HashMap<usize,Vec<usize>>, all_groups:&[String], uf:&UnionFind<usize>) -> HashMap<String, Vec<String>> {
    let mut merged_bparts:HashMap<String,Vec<String>> = HashMap::new();
    for (j, old_g) in all_groups.iter().enumerate(){
        let f_txp = old_g.clone().split("_").map(|x| x.parse::<usize>().unwrap()).collect::<Vec<usize>>()[0];
        let g_ind = uf.find(f_txp); //the vertex index
        let strings: Vec<String> = groups[&g_ind].iter().map(|n| n.to_string()).collect();
        let m_group = format!("{}", strings.join("_"));
        merged_bparts.entry(m_group).or_insert_with(Vec::new).push(all_groups[j].clone());
    }
    return merged_bparts
}

pub fn add_biparitions(g_bp_map:&HashMap<String, HashMap<String,u32>>, inds_bp_map:&HashMap<String, Vec<String>>) -> HashMap<String, HashMap<String, u32>> {
    let mut merged_bp_map:HashMap<String, HashMap<String, u32>> = HashMap::new();
    for (m_group, groups) in inds_bp_map.iter() {
        //println!("Merged {}", m_group);
        //println!("{:?}", groups);
        let m_bpart_key = merged_bp_map.entry(sort_group_id(&m_group.clone())).or_insert(HashMap::new());
        for g in groups.iter() {
            if !g_bp_map.contains_key(g){
                panic!("key {} not found", g);
            }
            for (b_part, count) in g_bp_map.get(g).unwrap().iter(){
                let c_count = m_bpart_key.entry(b_part.clone()).or_insert(0);
                *c_count += count;
            }
        }       
    }
    merged_bp_map
}