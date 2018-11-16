.. highlight:: shell

======================
Annotation Algorithm
======================

.. note:: There are several places where hard-coded logic was added to make the algorithm work with certain sequences. Hard coded logic is marked by **HARD CODED LOGIC** in the code. For instance, in :ref:`seq_search` I added logic to annotate exon 8 for class I last, because mapping it first causes issues due to the size of the feature.


#. **Check if locus is provided** (:ref:`blast`)

    * yes, then continue to step 2
    * no, then blast to get locus

#. **Check if any exact matches exist** (:ref:`ref`)

    * yes, then return annotation associated with the exact match and go to step 7
    * no, then go on to step 3

#. **Blast sequence and get list of alleles** (:ref:`blast`)

    * If all of the returned sequences are partial then the last sequence will be replaced with a fully characterized sequence.

#. **Iterate through the list and try to annotate with reference sequences** (:ref:`seq_search`)

    * Break reference up into features
    * Search for each feature in the provided sequence. When a feature is found record the coordinates and remove the mapped sequence from the unknown coordinates.
    * If all features are mapped, then go to step 7, else..
    * Try and assemble the remaining features based on what has already been mapped. Since we know the coordinates of the mapped features and the remaining unmapped sequence, we can determine if the unmapped sequences fall between two mapped features or at the ends/beginning after/before mapped sequences.
    * If all features could be mapped, then go to step 7 else go back to step 4A using any partial annotations for each reference sequence. If no annotation could be created after searching all of the reference sequences, then move on to step 5.
    * Partial annotations are :ref:`ann` objects with ``mapping``, ``blocks``, ``covered`` attributes. The ``mapping`` attribute is a dictionary with each position being a key and the values being features if they are mapped. The ``blocks`` attribute is a list of lists, which each list representing the positions of the unmapped parts of the sequence. If the whole sequence isn't mapped then there will only be one list (block).

#. **Loop through each reference sequence and do targeted alignments** (:ref:`align`)

    * Break up each reference sequence into features and create feature combos that will be used for doing alignments. Order the feature combos by the ones that make the most sense first.
    *  Do targeted alignments on all of the remaining blocks of sequences that have not yet been mapped.
  
        * If a high enough proportion of the unmapped sequence maps and the deletion/insertion rate is low enough, then extract the unmapped sequence from the alignment and map it.
        * If all features are mapped then go to step 7, else..
        * run step 4 with the updated partial annotation to see if the annotation can now be assembled. Go to step 7 if all features are mapped else..
    * Loop through all feature combinations for all reference sequences. This slows down the annotation if it's very novel. For instance, if it's a new feature sequence and that specific feature has only been reported in IMGT/HLA a few time for a given locus. The acceptance rate for the alignments is decreased slightly after each loop. For class I that decrease stops after the second reference sequence, but for class II it will keep going lower.
    * Rerun targeted alignment but with exons only combinations. 

#. **Do a full sequence alignment and use any partial annotation** (:ref:`align`)

    * If this fails and the rerun flag is set to ``True``, then rerun the whole annotation process starting from step 1. This time, skip the first reference allele that was used for doing the annotation and increase the number of reference alleles used by 1. 

#. **Generate GFE notation** (:ref:`gfe`)

    * Once a complete annotation is generated the GFE notation will be made
    * If the sequence only contains A,T,C or G, then a GFE notation can be created




















