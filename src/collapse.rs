use serde::{Serialize, Deserialize};
use petgraph::unionfind::UnionFind;
use std::collections::{HashMap, HashSet, BTreeMap};
use std::iter::FromIterator;
use crate::binary_tree::{sort_group_id};

fn create_union_find(g:&[String], ntxps:usize) -> UnionFind<usize> {
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

fn get_merged_bparts(groups:&HashMap<usize,Vec<usize>>, 
    all_groups_bpart:&HashMap<String, HashMap<String, u32>>, 
    uf:&UnionFind<usize>) -> HashMap<String, HashMap<String, u32>> {
    let all_groups:Vec<String> = all_groups_bpart.keys().cloned().collect();
    let mut merged_bparts:HashMap<String, HashMap<String,u32>> = HashMap::new();
    for (j, old_g) in all_groups.iter().enumerate(){
        let f_txp = old_g.clone().split("_").map(|x| x.parse::<usize>().unwrap()).collect::<Vec<usize>>()[0];
        let g_ind = uf.find(f_txp); //the vertex index
        let strings: Vec<String> = groups[&g_ind].iter().map(|n| n.to_string()).collect();
        let m_group = format!("{}", strings.join("_"));
        let m_bpart_key = merged_bparts.entry(sort_group_id(&m_group.clone())).or_insert(HashMap::new());

        for (b_part, count) in all_groups_bpart.get(&old_g.clone()).unwrap().iter(){
            let c_count = m_bpart_key.entry(b_part.clone()).or_insert(0);
            *c_count += count;
        }
    }
    return merged_bparts;
}

pub fn merge_groups(all_groups_bpart:&HashMap<String, HashMap<String, u32>>, 
                    ntxps:usize) -> HashMap<String, HashMap<String, u32>> {

    let all_groups:Vec<String> = all_groups_bpart.keys().cloned().collect();
    let g_union = create_union_find(&all_groups, ntxps as usize);
    let mut groups = HashMap::new();
    for i in 0..ntxps {
        let root = g_union.find(i);
        if root != i {
            groups.entry(root).or_insert(vec!(root)).push(i);
        }
    }
    get_merged_bparts(&groups, &all_groups_bpart, &g_union)
    
} 