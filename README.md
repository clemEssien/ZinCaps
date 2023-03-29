# ZinCaps

ZinCaps for Protein - zinc binding site predictionis implemented using Keras2.1.1 and Tensorflow. Installation has been tested in Linux and Mac OS X with Python 2.7 and 3.6 
# Training and testing data

Sample fasta files are provided for both training and testing i.e.
- sample-for-training.fasta
- sample-for-prediction.fasta

# Required Libraries
```
pip install pandas
pip install numpy
pip install scipy
pip install h5py==2.10.0
pip install -v keras==2.1.1
pip install tensorflow==1.15 (or GPU supported tensorflow, refer to https://www.tensorflow.org/install/ for instructions)
```
 - This is the Tensorflow version, you must change the backend to TensorFlow.
If you have run Keras at least once, you will find the Keras configuration file at:
$HOME/.keras/keras.json
If it isn’t there, you can create it. 
Change the default configuration file into:
```sh
{	
    "image_dim_ordering": "th",
    "epsilon": 1e-07,
    "floatx": "float32",
    "backend": "tensorflow"
}
```
# Running on GPU or CPU

>If you want to use GPU, you also need to install [CUDA]( https://developer.nvidia.com/cuda-toolkit) and [cuDNN](https://developer.nvidia.com/cudnn); refer to their websites for instructions. 
CPU is only suitable for prediction not training. 

For details of other parameters, run:
```sh
python train_models.py --help
```
or
```
python train_models.py -h
```

## Training
To train a custom model, run the following command:

``` python train_models.py -input [fasta file] -output-prefix [prefix of pre-trained model] -residue-types [candidate binding residues] - window 12 ```

* fasta file: A file in fasta format, a sample is provided
* residue-types: These refer to the candidate amino acid residues that bind to Zinc  
* window - refers to the number of amino acid residues to the left part or right part adjacent to the potential binding residue. An input sequence fragment used in this work was 25 i.e. (2 * (window size)+1). So the window size will be 12. 


### For example, 
``` python train_models.py -input sample-for-training.fasta -residue-types C,H,E,D -output-prefix models/output - window 12  ```

## Prediction
To run prediction, use the following command:
```sh
python predict.py -input [fasta file to run prediction on] -model-prefix [prefix of the path to the pre-trained model] -output [custom specified prefix for the prediction results] 
```

### For example,
``` python predict.py -input sample-for-prediction.fasta -residue-types C,H,E,D  -model-prefix models/output -output results ```

* The command above  would produce an output in a text file called results.txt
## Citation：
Please, if you use this code, cite the following paper:
C. Essien, D. Wang and D. Xu, "Capsule Network for Predicting Zinc Binding Sites in Metalloproteins," 2019 IEEE
International Conference on Bioinformatics and Biomedicine (BIBM), San Diego, CA, USA, 2019, pp. 2337-2341.

License
----
GNU v2.0
