use serde::{Serialize, Deserialize};
use std::collections::{HashMap, BTreeMap, HashSet};

fn print_type_of<T>(_: &T) {
    println!("{}", std::any::type_name::<T>());
}

#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct TreeNode {
    pub id: String,
    pub left: Option<Box<TreeNode>>,
    pub right: Option<Box<TreeNode>>,
}


impl TreeNode {
    pub fn create_leaf(id: String) -> TreeNode {
        TreeNode {
            id: id,
            left: None,
            right: None
        }
    }
    
    pub fn create_group(n1: TreeNode, n2: TreeNode) -> TreeNode {
        TreeNode {
            id:format!("{}_{}",n1.id, n2.id),
            left: Some(Box::new(n1)),
            right: Some(Box::new(n2))
        }
    }
    
    // Create group from t1_t2
    fn create_single_group(s:String) -> TreeNode {
        let txps: Vec<&str> = s.split("_").collect();
        let mut root = TreeNode::create_leaf(txps[0].to_string());
        
        for i in txps.iter().skip(1){
            let other = TreeNode::create_leaf(i.to_string());
            root = TreeNode::create_group(root, other); 
        }
        root
    }
    
    pub fn create_group_from_order(s:String) -> TreeNode {
        let txps: Vec<&str> = s.split("gr").collect();
        let mut root = TreeNode::create_single_group(txps[0].to_string());
        
        for i in txps.iter().skip(1) {
            let other = TreeNode::create_single_group(i.to_string());
            root = TreeNode::create_group(root, other); 
        }
        root
    }
    
    
    // Might Borrow the implementation defined in
    // https://sachanganesh.com/programming/graph-tree-traversals-in-rust/
    pub fn traverse_tree(&self) {
        if ! self.left.is_none(){
            println!("root is {}", self.id);
            let d = self.left.as_ref().unwrap();
            println!("left is {}", d.id);
            self.left.as_ref().unwrap().traverse_tree();
        }
        if ! self.right.is_none(){
            println!("root is {}", self.id);
            let d = self.right.as_ref().unwrap();
            println!("right is {}", d.id);
            self.right.as_ref().unwrap().traverse_tree();
        }
    }
}

// Takes the set of all txps in the tree, txps in the child node (edge from current root) and returns the corresponding bipartition as string
// Criteria - Each bipartition is first sorted and then returned as parent and child
pub fn get_bipart_split(par:&HashSet<u32>, child:&str)  -> String {
    let child_set: HashSet<u32> = child.split("_").map(|x| x.parse::<u32>().unwrap()).collect();
    //println!("{:?}",child_set);
    let mut req_par:Vec<u32> = par.difference(&child_set).cloned().collect();
    req_par.sort();
    //println!("{:?}", req_par);
    //https://stackoverflow.com/questions/53115999/what-is-the-idiomatic-way-of-converting-a-vec-of-references-to-a-vec-of-values
    //https://hermanradtke.com/2015/06/22/effectively-using-iterators-in-rust.html
    let mut child_set:Vec<u32> = child_set.iter().cloned().collect(); //Why iter not into_iter, clones
    child_set.sort();
    if req_par.len() == 1 && child_set.len() == 1 {
        req_par.push(child_set[0]);
        req_par.sort();
        return format!("{}bp{}", req_par[0].to_string(), req_par[1].to_string());
    }
    let (req_par, child_set) = if req_par.len() > child_set.len() {
        (req_par, child_set)
    } else {
        (child_set, req_par)
    };
    //vec![req_par, child_set]
    let par_str: String = req_par
                        .iter()
                        .map(|v| v.to_string())
                        .collect::<Vec<String>>()
                        .join("_");
    let child_str: String = child_set
                        .iter()
                        .map(|v| v.to_string())
                        .collect::<Vec<String>>()
                        .join("_");
    
    format!("{}bp{}", child_str, par_str)
}

// pub fn bipart_splits_writer(split_file: &mut File, split: String, ) -> Result<bool,io::Error> {
//     writeln!

// }

pub fn compute_bipart_count(node:&TreeNode, bp_map:&mut HashMap<String,u32>,
     dir_bp_map:&mut HashMap<String,u32>, root_set:&HashSet<u32>, g_bipart:&mut Vec<String>) {
    if ! node.left.is_none(){
        //println!("root is {}", node.id);
        let split = get_bipart_split(root_set, &node.left.as_ref().unwrap().id);
        if ! dir_bp_map.contains_key(&split) {
            g_bipart.push(split.clone());
            dir_bp_map.insert(split.clone(), 1);
            let count = bp_map.entry(split).or_insert(0);
            *count += 1;
        }
        //println!("left is {}", d.id);
        compute_bipart_count(node.left.as_ref().unwrap(), bp_map, dir_bp_map, root_set, g_bipart);
    }
    if ! node.right.is_none(){
        //println!("root is {}", node.id);
        let split = get_bipart_split(root_set, &node.right.as_ref().unwrap().id);
        if ! dir_bp_map.contains_key(&split) {
            g_bipart.push(split.clone());
            dir_bp_map.insert(split.clone(), 1);
            let count = bp_map.entry(split).or_insert(0);
            *count += 1;
        }
        compute_bipart_count(node.right.as_ref().unwrap(), bp_map, dir_bp_map, root_set, g_bipart);
    }
}

pub fn compute_bipart_count2(node:&TreeNode, bp_map:&mut HashMap<String,u32>,
    dir_bp_map:&mut HashMap<String,u32>) {
   if ! node.left.is_none(){
       //println!("root is {}", node.id);
       //let split = get_bipart_split(root_set, &node.left.as_ref().unwrap().id);
       let bpart = sort_group_id(&node.left.as_ref().unwrap().id);
       if dir_bp_map.contains_key(&bpart) {
           println!("bpart repeats in left");
       }
        //    g_bipart.push(split.clone());
           dir_bp_map.insert(bpart.clone(), 1);
           let count = bp_map.entry(bpart.clone()).or_insert(0);
           *count += 1;
       
       //println!("left is {}", d.id);
       //compute_bipart_count(node.left.as_ref().unwrap(), bp_map, dir_bp_map, root_set, g_bipart);
       compute_bipart_count2(node.left.as_ref().unwrap(), bp_map, dir_bp_map);
   }
   if ! node.right.is_none(){
       //println!("root is {}", node.id);
       //let split = get_bipart_split(root_set, &node.right.as_ref().unwrap().id);
       let bpart = sort_group_id(&node.right.as_ref().unwrap().id);
       if dir_bp_map.contains_key(&bpart) {
           println!("bpart repeats in right");  
       }
           //g_bipart.push(split.clone());
           dir_bp_map.insert(bpart.clone(), 1);
           let count = bp_map.entry(bpart.clone()).or_insert(0);
           *count += 1;
       
       //compute_bipart_count(node.right.as_ref().unwrap(), bp_map, dir_bp_map, root_set, g_bipart);
       compute_bipart_count2(node.right.as_ref().unwrap(), bp_map, dir_bp_map);
   }
}

pub fn sort_group_id(group: &str) -> String{
    let mut tps:Vec<u32> = group.split("_")
            .map(|x| x.parse::<u32>().unwrap())
            .collect();
    tps.sort();
    let mut group=tps.iter()
            .map(|x| x.to_string())
            .collect::<Vec<String>>()
            .join("_");
    group
}

pub fn get_binary_rooted_newick_string(node:&TreeNode) -> String {
    if node.left.is_none() && node.right.is_none() {
        return node.id.clone();
    }
    else {
        let l = get_binary_rooted_newick_string(node.left.as_ref().unwrap());
        let r = get_binary_rooted_newick_string(node.right.as_ref().unwrap());
        return format!("({},{})",l.clone(),r.clone());
    }
}
// fn main () {
//     let x = TreeNode {id:"12".to_string(), left: None, right: None};
//     // let y = TreeNode {id:"12".to_string(), left: None, right: None};
//     // // let mut z = TreeNode::create_group(x,y);
//     // println!("{}", z.id);
//     let s = String::from("6_7gr3_4gr1_2_5");
//     let d = TreeNode::create_group_from_order(s);
    
//     let par_set:HashSet<u32> = d.id.split("_").map(|x| x.parse::<u32>().unwrap()).collect();
    
//     // println!("{:?}",par_set);
//     // println!("{}",d.id);
//     // let e = TreeNode::get_bipart(&par_set, &d.left.unwrap().id);
//     // println!("{:?}", e);
//     //println!("{:?}", e[1]);
//     let mut bipart_counter: HashMap<String, u32> = HashMap::new();
//     compute_bipart_count(&d, &mut bipart_counter, &par_set);
    
//     //println!("{:?}", !x.left.is_none());
//     //let e = d.left.as_ref().unwrap();
//     // d.traverse_tree();
//     // d.traverse_tree();
//    // print_type_of(&d.left.as_ref().unwrap());
//     //if(x.left == None){
//         //println!("{}", d.left.unwrap().id);}
    
//     // let mut y = TreeNode {id:"123", left: None, right: None};
//     // let mut z = TreeNode {id:"123", left: Some(Box::new(x)), right: Some(Box::new(y))};
//     //let d = TreeNode::create_leaf("aa");
//     //let e = TreeNode::create_leaf("ee");
    
//     //println!("{}",f)
//     // x.insert("z");
//     // x.insert("b");
//     // x.insert("c");
//     // assert!(x == Node {
//     //     val: "m",
//     //     l: Some(Box::new(Node {
//     //         val: "b",
//     //         l: None,
//     //         r: Some(Box::new(Node { val: "c", l: None, r: None })),
//     //     })),
//     //     r: Some(Box::new(Node { val: "z", l: None, r: None })),
//     // });
// }

