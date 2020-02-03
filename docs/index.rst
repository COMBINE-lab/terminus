.. terminus documentation master file, created by
   sphinx-quickstart on Mon Feb  3 15:38:59 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to terminus's documentation!
====================================

What is Terminus
----------------
Terminus is a tool that implements a new method for analyzing transcript-level 
abundance estimates from RNA-seq data. Terminus works downstream of salmon, and 
collapses individual transcripts into groups whose total transcriptional output 
can be estimated more accurately and robustly.

Requirements
------------
Terminus uses the [cargo](https://github.com/rust-lang/cargo) build system and package manager.
To build terminus from source, you will need to have rust (ideally v1.40 or greater) installed. 
Then, you can build terminus by executing:

.. code-block:: bash

    cargo build --release

Input to terminus
-----------------
Terminus expects salmon to be run on the raw fastq files with two non-default option enabled, 
``--numGibbsSamples`` and ``-d``. 
A typical run of salmon that is ideal for terminus is as follows,

.. code-block:: bash

    salmon quant -la -i <index> -1 <fatsq_file(s)> -2 <fatsq_file(s)> --numGibbsSamples 100 -d -o <output>

This step can be run on multiple samples in case, one wish to run terminus on multiple samples.

Terminus Grouping
-----------------
At the first step terminus groups the transcripts per experiment. This step deals with each experiment
independently. To run grouping terminus has to be run with the following command,

.. code-block:: bash

    target/release/terminus group -m <> --tolerance <> -d input_dir/<experiment> -o output_dir

The above command expects ``input_dir/<experiment>`` to contain all the files that are written by salmon
along with the ``bootstraps`` directory and the equivalence class file. Note that it *assumes* the ``experiment``
is the experiment name inside the `input_dir`. Terminus would create a direcoty named ``experiment`` in
``<output_dir>``. After a successfull run there will be a ``groups.txt`` file written in ``output_dir/experiment``.

This step is parallalizable trivially by taking help of ``gnu parallel`` or such tools.

Grouping options
----------------
* ``-d``: Input directory where the salmon files are written.
* ``-m``: The threshold value for passing the `min-spread`, for a posterior distribution, spread is defined as (max - min)/mean value. 
* ``-o``: output directory where the a folder by the root name of experiment would be created.
* ``--tolerance``: The allowable difference between the weight vectors for transcripts to consider them as identical.
* ``--seed``: An iteger value for deciding the seed. This seed would be used for all the random number generations. The default seed is `10`. 


Terminus Collapsing
-------------------
Given the groups are written in the directory as specified above, collapsing operation collapses the groups
and create a `consensus` group that is common across all the experiments.

.. code-block:: bash

    target/release/terminus collapse -c <>  -d <input_dir/experiment_1> <input_dir/experument_2> ... -o output_dir

The results of this step would be written in `output_dir` in the corresponding folders `experiment_1`, `experument_2`
etc.

Collapsing options
------------------
* ``-d``: Input directory where all the root level directories of salmon is present
* ``-c``: The consensus threshold, determines the number of experiments a group has to appear into. A value of `0.5` dictates, that the final grups are at least present in half of the experiments.


Output
------
At the end of the above two steps the final directory would contain a the experiments directory with
in the *exact* same way salmon writes output. Although the number of transcripts in the final output
would be `total_number_of_transcripts - transcripts_collapsed + groups`. 


Contents:

.. toctree::
   :maxdepth: 2



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
