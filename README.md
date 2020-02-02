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

Terminus has two sub-commands, `group` and `collapse`.
