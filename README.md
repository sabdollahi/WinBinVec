<h1> WinBinVec </h1>
In this repository, we implemented seven different deep learning models for cancer type prediction and metastasis prediction using protein-protein interactions' sequences:
<ul>
  <li>WinBinVec (Window-based Binary Vectors)</li> 
  <li>WinSelfAtt (Window-based Self Attention)</li> 
  <li>BindingAffinity (Binding affinity changes upon mutations)</li>
  <li>gGapDipeptide (g gap dipeptide modeling of PPIs)</li> 
  <li>SMFM (Single Matrix Factorization Model)</li>
  <li>ExprGCNPPI (Gene Expression-based Graph Convolution Network)</li>
  <li>PhyChemPPI (Physicochemical Properties of PPIs)</li>
</ul>
Besides, we designed a model, <em>FusionPPI</em>, that utilizes the outperformed feature sets, namely window-based binary vectors, degree of pathogenicity, binding affinity changes upon mutations, and gene expressions. <br>
You can download the gene expression data from the following link:
https://drive.google.com/file/d/11pslP-6fsuj8C-9vSY2cOp3QYkcz2b5i/view?usp=sharing
Furthermore, this repository provides the codes for data processing and data cleaning including: PDB-related processing, Standard deviation of the mutation accumulation, Extraction of the positions of the available amino acids in each PDB files, Window selection method, and plotting the figures.
