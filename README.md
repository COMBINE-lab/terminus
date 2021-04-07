![CI](https://github.com/COMBINE-lab/Terminus/workflows/CI/badge.svg?branch=master)

<p align="center">
[<img src="https://upload.wikimedia.org/wikipedia/commons/thumb/d/dd/Design_for_a_Stained_Glass_Window_with_Terminus%2C_by_Hans_Holbein_the_Younger.jpg/800px-Design_for_a_Stained_Glass_Window_with_Terminus%2C_by_Hans_Holbein_the_Younger.jpg" alt="drawing" width="200">]
</p>

What is terminus?
=================

Terminus is a program for analyzing transcript-level abundance estimates from RNA-seq data, 
computed using [salmon](https://github.com/COMBINE-lab/salmon), and collapsing individual transcripts 
into groups whose total transcriptional output can be estimated more accurately and robustly.

The groups computed by terminus represent abundance estimation reported at the resolution that 
is actually supported by the underlying experimental data. In a typical experiment, this is 
neither at the gene level nor the transcript level. Some transcripts, even from complex, multi-isoform genes, 
can have their abundances confidently estimated, while other transcripts cannot. Rather than pre-defining 
the resolution at which the analysis will be performed, and subjecting the results to unnecessary uncertainty 
or insufficient biological resolution, terminus allows the determination of transcriptional groups that can 
be confidently ascertained in a given sample, and represents, in this sense, a data-driven approach to 
transcriptome analysis.


How to build terminus
---------------------

Terminus uses the [cargo](https://github.com/rust-lang/cargo) build system and package manager.  To build terminus from source, you will need to have rust (ideally v1.40 or greater) installed.  Then, you can build terminus by executing:

```
$ cargo build --release
```

from the top-level directory of this repository.  This will produce an executable in `target/release/terminus`.

How to use terminus
-------------------

Terminus has two sub-commands, `group` and `collapse`. For detailed tutorial and usage please visit the [tutorial](https://combine-lab.github.io/terminus-tutorial/2020/running-terminus/) or the [documentation](https://terminus.readthedocs.io/en/latest/).

Updates
-------
`group` and `collapse` step are now updated where `group` creates a tree for each group and `collapse` creates either a consensus tree directly or the input that could be fed to a consenus tree algorithm. Under `collapse` there are additional sub-arguments - `m` which tells whether the overlapping groups across samples are to be merged (deafult true) and `merge_type` (default phylip) indicates the type of consensus method - majority consensus or the majority rule extended method. For more details for these refer to [this](https://evolution.genetics.washington.edu/phylip/doc/consense.html)

The majority consensus is supported by biopy package in python2. Our approach would produce the file ```cluster_bipart_splits.txt``` which would be used as input for biopy. This is obtained by using `BP` as the argument to the `merge_type` parameter. The majority rule extended method is being computed via phylip and the collapse step would also yield the final consensus tree in `cluster_nwk.txt`. A limitation with running phylip option is that currently it is only implemented for the flag `true` for the parameter `m` unlike `BP` which is independent of it.

An example of code snippet
```
$target/release/terminus collapse -c 0  -d data/Salmon_shuffled/Sample1 data/Salmon_shuffled/Sample2 data/Salmon_shuffled/Sample3 data/Salmon_shuffled/Sample4 data/Salmon_shuffled/Sample5 data/Salmon_shuffled/Sample6 -o data/term -m true --merge_type phylip
```

Things in the process
<ol>
  <li> Building majority rule consensus tree using Phylip. </li>
  <li> Support in phylip which can take both merged and not merged groups as inputs. </li>
  <li> Replacing the index names in the trees with the actual transcript names for Phylip.</li>
  <li> Calling the C function for consensus directly from rust. The current process involves calling phylip consensus algorithm passing an input file for group containing all trees which invovles a lot of IOs that makes the process slow.</li>
</ol>

Authors
-------

Hirak Sarkar, Avi Srivastava, H&egrave;ctor Corrada Bravo, Michael I. Love, Rob Patro
