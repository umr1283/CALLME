
The example data in example_FDT.txt is composed of 12 samples, of which six have clonal mosaic events.

NOTE : as cluster centers (mean value of R and theta for each rs and each genotype) are only calculated on 12 samples,
the normalization step is not perfect. But it still works as an example.

Commands to process from here:


	 > cd ..
	 > perl CALLME.pl --data_file=example/example_FDT.txt --n_markers=80000 --rs_only=true --n_procs=2 --T=2 --aAlpha=0.1 --min_seg_length=75

Here, at each step of analysis, two processes will be launched in parallel. You can adjust this paramater according to you system.
