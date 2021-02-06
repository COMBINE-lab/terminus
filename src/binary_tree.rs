fn print_type_of<T>(_: &T) {
    println!("{}", std::any::type_name::<T>());
}

struct TreeNode {
    id: String,
    left: Option<Box<TreeNode>>,
    right: Option<Box<TreeNode>>,
}


impl TreeNode {
    fn create_leaf(id: String) -> TreeNode {
        TreeNode {
            id: id,
            left: None,
            right: None
        }
    }
    
    fn create_group(n1: TreeNode, n2: TreeNode) -> TreeNode {
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
        let txps: Vec<&str> = s.rsplit("gr").collect();
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
fn main () {
    let x = TreeNode {id:"12".to_string(), left: None, right: None};
    // let y = TreeNode {id:"12".to_string(), left: None, right: None};
    // // let mut z = TreeNode::create_group(x,y);
    // println!("{}", z.id);
    let s = "6_7gr3_4gr1_2_5".to_string();
    let d = TreeNode::create_group_from_order(s);
    //println!("{:?}", !x.left.is_none());
    //let e = d.left.as_ref().unwrap();
    d.traverse_tree();
    d.traverse_tree();
   // print_type_of(&d.left.as_ref().unwrap());
    //if(x.left == None){
        //println!("{}", d.left.unwrap().id);}
    
    // let mut y = TreeNode {id:"123", left: None, right: None};
    // let mut z = TreeNode {id:"123", left: Some(Box::new(x)), right: Some(Box::new(y))};
    //let d = TreeNode::create_leaf("aa");
    //let e = TreeNode::create_leaf("ee");
    
    //println!("{}",f)
    // x.insert("z");
    // x.insert("b");
    // x.insert("c");
    // assert!(x == Node {
    //     val: "m",
    //     l: Some(Box::new(Node {
    //         val: "b",
    //         l: None,
    //         r: Some(Box::new(Node { val: "c", l: None, r: None })),
    //     })),
    //     r: Some(Box::new(Node { val: "z", l: None, r: None })),
    // });
}
