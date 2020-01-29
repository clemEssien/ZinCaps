# ActiveSitePrediction

This code was designed to run on GPU 
 
Instructions
1. Run /ZinCaps/ActiveSitePrediction/src/DataProcessing/genSequence.py
to generate annotated protein sequences and also to generate 5 fold cross validation.

2. To train, run the following command

    python train_models.py -input "/ZinCaps/ActiveSitePrediction/lib/K-Fold/annotated_sequence.fasta_training_annotated_4.fasta" 
    -output-prefix /ZinCaps/ActiveSitePrediction/data/weights/ZincCaps_4_model -checkpointweights 
    /home/ucexbw/ZinCaps/ActiveSitePrediction/data/weights/ZincCaps_4_weights -residue-types C,H,E,D -nclass=5 -maxneg 30 -window 25

3. To run prediction, run the following
    python predict_and_evaluate_sharedarch.py -residue-types C,H,E,D -window 25 -model-prefix 
    /ZinCaps/ActiveSitePrediction/data/weights/ZinCaps_ 

4. The output would be in a textfile located in /ZinCaps/ActiveSitePrediction/data/output/

5. To compare performance with CNN, code used is available at 
     https://github.com/duolinwang/MusiteDeep/tree/master/MusiteDeep_Keras2.0
