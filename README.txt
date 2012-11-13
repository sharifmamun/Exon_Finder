A novel approach to find the intron-exon positions in EST contigs of Fly insects. 
This codes will need specific datasets to run. If you need datasets, please email me.

The sequences of the code execution will be:
  1. getesty.py
  2. result_handler.py
  3. frame_mapper.py

It's possible to run all the scripts from the same file but for safety and dependency with one another, I myself run them one after another.

* checker.py files should exist in the same directory!!!

Change needed!!! (Till everything is automated...)
result_hnadler.py
        | In main(): change the name of handle and hadle1 file names
        | temp = grep("Ceraphronidae_sp", read_data)
        | "Ceraphronidae_sp" should be replaced by the related string of the species.

Change in helping script:
-------------------------
checker.check.py
        | make sure you are looking for the perfect file: in handle2 variable!