<h1> WinBinVec </h1>
In this repository, we implemented seven different deep learning models for cancer type prediction and metastasis prediction using protein-protein interactions' sequences:
<ul>
  <li><em>WinBinVec</em> (Window-based Binary Vectors)</li> 
  <li><em>WinSelfAtt</em> (Window-based Self Attention)</li> 
  <li><em>BindingAffinity</em> (Binding affinity changes upon mutations)</li>
  <li><em>gGapDipeptide</em> (g gap dipeptide modeling of PPIs)</li> 
  <li><em>SMFM</em> (Single Matrix Factorization Model)</li>
  <li><em>ExprGCNPPI</em> (Gene Expression-based Graph Convolution Network)</li>
  <li><em>PhyChemPPI</em> (Physicochemical Properties of PPIs)</li>
</ul>
Besides, we designed a model, <em>FusionPPI</em>, that utilizes the outperformed feature sets, namely window-based binary vectors, degree of pathogenicity, binding affinity changes upon mutations, and gene expressions. <br>
To run each of the above models, use the following instructions: <br>
<ol>
  <li> Download all the input feature sets and the necessary files. </li> 
  <li> Create a directory with name "DATASET" beside the models' files. </li>
  <li> Extract the downloaded files into the "DATASET" directory. </li>
</ol>  
<br>
You can download the necessary datasets from the following links:
<br>
<b> WinBinVec Input Features: </b>
https://drive.google.com/file/d/1eALbG7wyZwjGsLG_XCy8QqtL40KkQiSp/view?usp=sharing
<br>
<b> Gene Expression Data: </b>
https://drive.google.com/file/d/11pslP-6fsuj8C-9vSY2cOp3QYkcz2b5i/view?usp=sharing
<br>
<b> g-Gap Dipeptide Input Features: </b>
https://drive.google.com/file/d/1Bo9x_dukeCviuhf9M-h7Oqk9YAMIt1gV/view?usp=sharing
<br>
<b> Physicochemical Properties Input Features: </b>
https://drive.google.com/file/d/1UHLZrPDPrfUS69OySL98af-_-mJN6_iP/view?usp=sharing
<br>
<b> WinSelfAtt Input Features: </b>
https://drive.google.com/file/d/1jNTfaVwXpichkNRSSscoktdlpmCcONQT/view?usp=sharing
<br>
<b> The other Input Features: </b>
https://drive.google.com/file/d/1vtrJL84HhCzqMAstyDW2_h5ffELptW8i/view?usp=sharing
<br>




Furthermore, this repository provides the codes for data processing and data cleaning including: PDB-related processing, Standard deviation of the mutation accumulation, Extraction of the positions of the available amino acids in each PDB files, Window selection method, and plotting the figures.
