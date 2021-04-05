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
<h2> Prerequisite </h2>
The following libraries are prerequisite for running the WinBinVec models:
<ul>
  <li> PyTorch &#8805; 1.7 </li>
  <li> pickle &#8805; 4.0 </li>
  <li> pandas &#8805; 1.1.5 </li>
  <li> sklearn &#8805; 0.19.1 </li>
  <li> numpy &#8805; 1.19.5 </li>
  <li> seaborn &#8805; 0.9.0 </li>
  <li> json &#8805; 2.0.9 </li>
  <li> matplotlib &#8805; 3.3.4 </li>
  <li> lifelines &#8805; 0.24.9 </li>
</ul>
<h2> Run the models </h2>
To run each of the above models, use the following instructions: <br>
<ol>
  <li> Download all the input feature sets and the necessary files. </li> 
  <li> Create a directory with name "DATASET" beside the models' files. </li>
  <li> Extract the downloaded files into the "DATASET" directory. </li>
    <li> Run: python3 WinBinVec.py</li>
</ol>  
If you want to get the other models' performance, you need to run <p  style="color:red;">python3 [A-MODEL-NAME].py</p> by replacing [A-MODEL-NAME] with the name of your model of interest.
<br>
We also designed a second version of the WinBinVec model, WinBinVec-FC. In this case, we utilized fully-connected layers modules instead of 1-dimensional convolutional module. The results (accuracy, PR, and AUC) reveal that the original version of the WinBinVec model (using 1-dimensional convolutional module) outperforms the second version.
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
<b> The other Input Features and the necessary files: </b>
https://drive.google.com/file/d/1vtrJL84HhCzqMAstyDW2_h5ffELptW8i/view?usp=sharing
<br>
<b> The preprocessing necessary files: </b>
https://drive.google.com/file/d/1iWQIPi02xnZ5ayU6f8kZvqobSiGjTJW5/view?usp=sharing
<br>



Furthermore, this repository provides the codes for data processing and data cleaning including: PDB-related processing, Standard deviation of the mutation accumulation, Extraction of the positions of the available amino acids in each PDB files, Window selection method, and plotting the figures. <br>
You can obtain the predicted essentiality values (cancer-associated PPIs) of the WinBinVec model and the MEDICI model using the following link: <br>
https://drive.google.com/file/d/1XSbqIRLFOXEpOHRUlxXk_OGHjxiM4o-0/view?usp=sharing
<br>

You can obtain the protein-protein interactions list including their PDB ids in the following box:
<table style="width:100%">
  <tr>
    <th>PPI Number</th>
    <th>Protein RefSeq Id 1</th>
    <th>PDB Id and Chain 1</th> 
    <th>Protein RefSeq Id 2</th>
    <th>PDB Id and Chain 2</th>
  </tr>
  <tr><td>1</td><td>NM_000020</td><td>4FAO_C</td><td>NM_001106</td><td>4FAO_E</td></tr>
<tr><td>2</td><td>NM_004757</td><td>4R3Z_A</td><td>NM_002887</td><td>4R3Z_B</td></tr>
<tr><td>3</td><td>NM_005163</td><td>3MV5_A</td><td>GSK3B</td><td>3MV5_C</td></tr>
<tr><td>4</td><td>NM_001626</td><td>1O6K_A</td><td>GSK3B</td><td>1O6K_C</td></tr>
<tr><td>5</td><td>NM_000477</td><td>4K71_A</td><td>NM_004107</td><td>4K71_B</td></tr>
<tr><td>6</td><td>NM_000477</td><td>4K71_A</td><td>NM_001136019</td><td>4K71_B</td></tr>
<tr><td>7</td><td>NM_000477</td><td>4N0U_D</td><td>NM_004048</td><td>4N0U_B</td></tr>
<tr><td>8</td><td>NM_013367</td><td>5KHU_I</td><td>NM_003903</td><td>5KHU_J</td></tr>
<tr><td>9</td><td>NM_173473</td><td>5KHU_E</td><td>NM_001256</td><td>5KHU_F</td></tr>
<tr><td>10</td><td>NM_000039</td><td>2MSE_A</td><td>NM_001654</td><td>2MSE_D</td></tr>
<tr><td>11</td><td>NM_000039</td><td>2MSE_A</td><td>NM_001256197</td><td>2MSE_D</td></tr>
<tr><td>12</td><td>NM_130384</td><td>5YZ0_C</td><td>NM_001184</td><td>5YZ0_A</td></tr>
<tr><td>13</td><td>NM_000489</td><td>5GRQ_C</td><td>NM_001141969</td><td>5GRQ_A</td></tr>
<tr><td>14</td><td>NM_000489</td><td>5GRQ_C</td><td>NM_001350</td><td>5GRQ_A</td></tr>
<tr><td>15</td><td>NM_000332</td><td>4J2L_A</td><td>NM_015125</td><td>4J2L_C</td></tr>
<tr><td>16</td><td>NM_003600</td><td>5ODT_A</td><td>NM_006342</td><td>5ODT_B</td></tr>
<tr><td>17</td><td>NM_021913</td><td>4RA0_C</td><td>NM_000820</td><td>4RA0_A</td></tr>
<tr><td>18</td><td>NM_014417</td><td>6QFM_B</td><td>NM_021960</td><td>6QFM_A</td></tr>
<tr><td>19</td><td>NM_005872</td><td>5MQF_K</td><td>NM_002669</td><td>5MQF_D</td></tr>
<tr><td>20</td><td>NM_005872</td><td>5MQF_K</td><td>NM_001253</td><td>5MQF_L</td></tr>
<tr><td>21</td><td>NM_000709</td><td>1V11_A</td><td>NM_000056</td><td>1V11_B</td></tr>
<tr><td>22</td><td>NM_004326</td><td>2VPB_B</td><td>NM_015617</td><td>2VPB_A</td></tr>
<tr><td>23</td><td>NM_004326</td><td>2GL7_C</td><td>NM_030756</td><td>2GL7_B</td></tr>
<tr><td>24</td><td>NM_138578</td><td>3R85_A</td><td>NM_014320</td><td>3R85_E</td></tr>
<tr><td>25</td><td>NM_138578</td><td>4HNJ_A</td><td>NM_014417</td><td>4HNJ_C</td></tr>
<tr><td>26</td><td>NM_017745</td><td>6B7G_B</td><td>NM_004529</td><td>6B7G_A</td></tr>
<tr><td>27</td><td>NM_017745</td><td>4HPL_A</td><td>NM_032673</td><td>4HPL_B</td></tr>
<tr><td>28</td><td>NM_021946</td><td>5JH5_D</td><td>NM_032590</td><td>5JH5_A</td></tr>
<tr><td>29</td><td>NM_197966</td><td>2KBW_B</td><td>NM_021960</td><td>2KBW_A</td></tr>
<tr><td>30</td><td>NM_001166</td><td>6HPR_A</td><td>NM_181838</td><td>6HPR_C</td></tr>
<tr><td>31</td><td>NM_001165</td><td>3M0A_D</td><td>NM_021138</td><td>3M0A_A</td></tr>
<tr><td>32</td><td>NM_139317</td><td>4AUQ_B</td><td>NM_021009</td><td>4AUQ_C</td></tr>
<tr><td>33</td><td>NM_004329</td><td>3QB4_B</td><td>GDF5</td><td>3QB4_A</td></tr>
<tr><td>34</td><td>NM_004329</td><td>2QJ9_D</td><td>NM_001200</td><td>2QJ9_B</td></tr>
<tr><td>35</td><td>NM_004329</td><td>2H62_C</td><td>NM_001106</td><td>2H62_D</td></tr>
<tr><td>36</td><td>NM_007294</td><td>1JM7_A</td><td>NM_000465</td><td>1JM7_B</td></tr>
<tr><td>37</td><td>NM_000059</td><td>1N0W_B</td><td>NM_002875</td><td>1N0W_A</td></tr>
<tr><td>38</td><td>NM_000059</td><td>1N0W_B</td><td>NM_133487</td><td>1N0W_A</td></tr>
<tr><td>39</td><td>NM_001211</td><td>5KHU_Q</td><td>NM_002358</td><td>5KHU_T</td></tr>
<tr><td>40</td><td>NM_003910</td><td>5MQF_Q</td><td>NM_016403</td><td>5MQF_R</td></tr>
<tr><td>41</td><td>NM_004048</td><td>1GZP_B</td><td>NM_001764</td><td>1GZP_A</td></tr>
<tr><td>42</td><td>NM_004048</td><td>1XZ0_B</td><td>NM_001763</td><td>1XZ0_A</td></tr>
<tr><td>43</td><td>NM_001743</td><td>2X0G_B</td><td>NM_004938</td><td>2X0G_A</td></tr>
<tr><td>44</td><td>NM_001228</td><td>3H11_B</td><td>NM_003879</td><td>3H11_A</td></tr>
<tr><td>45</td><td>NM_032996</td><td>1NW9_B</td><td>NM_001167</td><td>1NW9_A</td></tr>
<tr><td>46</td><td>NM_053056</td><td>2W96_A</td><td>NM_000075</td><td>2W96_B</td></tr>
<tr><td>47</td><td>NM_001238</td><td>5L2W_B</td><td>NM_001798</td><td>5L2W_A</td></tr>
<tr><td>48</td><td>NM_001256</td><td>5KHU_F</td><td>NM_004661</td><td>5KHU_C</td></tr>
<tr><td>49</td><td>NM_001791</td><td>1GRN_A</td><td>NM_004308</td><td>1GRN_B</td></tr>
<tr><td>50</td><td>NM_001254</td><td>4I5L_B</td><td>NM_014225</td><td>4I5L_A</td></tr>
<tr><td>51</td><td>NM_004360</td><td>3FF7_A</td><td>NM_005810</td><td>3FF7_C</td></tr>
<tr><td>52</td><td>NM_004064</td><td>1JSU_C</td><td>NM_001798</td><td>1JSU_A</td></tr>
<tr><td>53</td><td>NM_004064</td><td>6ATH_C</td><td>NM_001237</td><td>6ATH_B</td></tr>
<tr><td>54</td><td>NM_000077</td><td>1BI7_B</td><td>CDK6</td><td>1BI7_A</td></tr>
<tr><td>55</td><td>NM_016507</td><td>4CXA_A</td><td>NM_001099402</td><td>4CXA_B</td></tr>
<tr><td>56</td><td>NM_001798</td><td>1FQ1_B</td><td>NM_005192</td><td>1FQ1_A</td></tr>
<tr><td>57</td><td>NM_001798</td><td>1FQ1_B</td><td>NM_001330173</td><td>1FQ1_A</td></tr>
<tr><td>58</td><td>NM_001798</td><td>1JSU_A</td><td>NM_001237</td><td>1JSU_B</td></tr>
<tr><td>59</td><td>NM_000075</td><td>3G33_A</td><td>NM_001760</td><td>3G33_B</td></tr>
<tr><td>60</td><td>NM_000075</td><td>5FWK_K</td><td>NM_007065</td><td>5FWK_E</td></tr>
<tr><td>61</td><td>CDK6</td><td>1BI8_A</td><td>NM_079421</td><td>1BI8_B</td></tr>
<tr><td>62</td><td>NM_014143</td><td>5IUS_C</td><td>NM_005018</td><td>5IUS_A</td></tr>
<tr><td>63</td><td>NM_000735</td><td>1XWD_A</td><td>NM_000510</td><td>1XWD_B</td></tr>
<tr><td>64</td><td>NM_004236</td><td>4D10_B</td><td>NM_016129</td><td>4D10_D</td></tr>
<tr><td>65</td><td>NM_006837</td><td>4D10_E</td><td>NM_006833</td><td>4D10_F</td></tr>
<tr><td>66</td><td>NM_006710</td><td>4D10_H</td><td>NM_001164095</td><td>4D10_G</td></tr>
<tr><td>67</td><td>NM_006710</td><td>4D10_H</td><td>NM_001164093</td><td>4D10_G</td></tr>
<tr><td>68</td><td>NM_006710</td><td>4D10_H</td><td>NM_016319</td><td>4D10_G</td></tr>
<tr><td>69</td><td>NM_006710</td><td>4D10_H</td><td>NM_001164094</td><td>4D10_G</td></tr>
<tr><td>70</td><td>NM_005211</td><td>4WRL_A</td><td>NM_000757</td><td>4WRL_B</td></tr>
<tr><td>71</td><td>NM_005211</td><td>4DKD_C</td><td>NM_152456</td><td>4DKD_A</td></tr>
<tr><td>72</td><td>NM_001895</td><td>1JWH_A</td><td>NM_001320</td><td>1JWH_C</td></tr>
<tr><td>73</td><td>NM_001098210</td><td>1JDH_A</td><td>NM_030756</td><td>1JDH_B</td></tr>
<tr><td>74</td><td>NM_001098210</td><td>3TX7_A</td><td>NM_205860</td><td>3TX7_B</td></tr>
<tr><td>75</td><td>NM_001098210</td><td>2GL7_A</td><td>NM_004326</td><td>2GL7_C</td></tr>
<tr><td>76</td><td>NM_001098210</td><td>1T08_A</td><td>NM_020248</td><td>1T08_B</td></tr>
<tr><td>77</td><td>NM_003592</td><td>1U6G_A</td><td>NM_014248</td><td>1U6G_B</td></tr>
<tr><td>78</td><td>NM_020640</td><td>3TDU_A</td><td>NM_003592</td><td>3TDU_C</td></tr>
<tr><td>79</td><td>NM_001291628</td><td>5ZAL_A</td><td>TARBP2</td><td>5ZAL_B</td></tr>
<tr><td>80</td><td>NM_177438</td><td>5ZAL_A</td><td>TARBP2</td><td>5ZAL_B</td></tr>
<tr><td>81</td><td>NM_030621</td><td>5ZAL_A</td><td>TARBP2</td><td>5ZAL_B</td></tr>
<tr><td>82</td><td>NM_001271282</td><td>5ZAL_A</td><td>TARBP2</td><td>5ZAL_B</td></tr>
<tr><td>83</td><td>NM_022552</td><td>5YX2_A</td><td>NM_013369</td><td>5YX2_B</td></tr>
<tr><td>84</td><td>NM_005228</td><td>5WB7_A</td><td>NM_001432</td><td>5WB7_E</td></tr>
<tr><td>85</td><td>NM_005228</td><td>3NJP_A</td><td>EGF</td><td>3NJP_C</td></tr>
<tr><td>86</td><td>NM_005228</td><td>1MOX_A</td><td>NM_003236</td><td>1MOX_C</td></tr>
<tr><td>87</td><td>NM_005228</td><td>4ZJV_A</td><td>NM_018948</td><td>4ZJV_C</td></tr>
<tr><td>88</td><td>NM_001430</td><td>6CZW_A</td><td>NM_001668</td><td>6CZW_B</td></tr>
<tr><td>89</td><td>NM_004431</td><td>3MBW_A</td><td>NM_004428</td><td>3MBW_B</td></tr>
<tr><td>90</td><td>NM_001982</td><td>4RIW_A</td><td>NM_005228</td><td>4RIW_B</td></tr>
<tr><td>91</td><td>NM_001042599</td><td>3U7U_A</td><td>NM_004495</td><td>3U7U_G</td></tr>
<tr><td>92</td><td>NM_001042599</td><td>3U7U_A</td><td>NM_013964</td><td>3U7U_G</td></tr>
<tr><td>93</td><td>NM_000122</td><td>5IVW_V</td><td>NM_000400</td><td>5IVW_W</td></tr>
<tr><td>94</td><td>NM_000082</td><td>4A11_B</td><td>NM_001923</td><td>4A11_A</td></tr>
<tr><td>95</td><td>NM_004456</td><td>5WUK_B</td><td>NM_003797</td><td>5WUK_A</td></tr>
<tr><td>96</td><td>NM_004456</td><td>6C23_K</td><td>NM_004973</td><td>6C23_E</td></tr>
<tr><td>97</td><td>NM_001257069</td><td>5V4B_B</td><td>NM_170679</td><td>5V4B_A</td></tr>
<tr><td>98</td><td>NM_033632</td><td>5V4B_B</td><td>NM_170679</td><td>5V4B_A</td></tr>
<tr><td>99</td><td>NM_001349798</td><td>5V4B_B</td><td>NM_170679</td><td>5V4B_A</td></tr>
<tr><td>100</td><td>NM_023110</td><td>5W59_B</td><td>NM_002010</td><td>5W59_A</td></tr>
<tr><td>101</td><td>NM_023110</td><td>1CVS_C</td><td>NM_002006</td><td>1CVS_A</td></tr>
<tr><td>102</td><td>NM_023110</td><td>1EVT_C</td><td>NM_000800</td><td>1EVT_A</td></tr>
<tr><td>103</td><td>NM_000141</td><td>2FDB_P</td><td>NM_006119</td><td>2FDB_M</td></tr>
<tr><td>104</td><td>NM_000141</td><td>2FDB_P</td><td>NM_033165</td><td>2FDB_M</td></tr>
<tr><td>105</td><td>NM_000141</td><td>4J23_A</td><td>NM_000800</td><td>4J23_B</td></tr>
<tr><td>106</td><td>NM_005117</td><td>6NFJ_C</td><td>NM_175737</td><td>6NFJ_A</td></tr>
<tr><td>107</td><td>NM_000141</td><td>1NUN_B</td><td>NM_004465</td><td>1NUN_A</td></tr>
<tr><td>108</td><td>NM_000142</td><td>1RY7_B</td><td>NM_000800</td><td>1RY7_A</td></tr>
<tr><td>109</td><td>NM_002006</td><td>1II4_A</td><td>NM_000141</td><td>1II4_E</td></tr>
<tr><td>110</td><td>NM_001110556</td><td>2MTP_A</td><td>NM_000419</td><td>2MTP_B</td></tr>
<tr><td>111</td><td>NM_002019</td><td>5T89_X</td><td>NM_001025366</td><td>5T89_V</td></tr>
<tr><td>112</td><td>FLT3</td><td>3QS7_E</td><td>NM_001204503</td><td>3QS7_A</td></tr>
<tr><td>113</td><td>FLT3</td><td>3QS7_E</td><td>NM_001204502</td><td>3QS7_A</td></tr>
<tr><td>114</td><td>FLT3</td><td>3QS7_E</td><td>NM_001459</td><td>3QS7_A</td></tr>
<tr><td>115</td><td>NM_002020</td><td>4BSK_A</td><td>NM_005429</td><td>4BSK_C</td></tr>
<tr><td>116</td><td>NM_182925</td><td>4BSK_A</td><td>NM_005429</td><td>4BSK_C</td></tr>
<tr><td>117</td><td>NM_001354989</td><td>4BSK_A</td><td>NM_005429</td><td>4BSK_C</td></tr>
<tr><td>118</td><td>NM_052905</td><td>4YC7_B</td><td>NM_001791</td><td>4YC7_A</td></tr>
<tr><td>119</td><td>NM_212476</td><td>2OCF_D</td><td>NM_000125</td><td>2OCF_A</td></tr>
<tr><td>120</td><td>NM_212476</td><td>5DC0_A</td><td>NM_007313</td><td>5DC0_B</td></tr>
<tr><td>121</td><td>NM_002015</td><td>4LG0_A</td><td>NM_005238</td><td>4LG0_B</td></tr>
<tr><td>122</td><td>NM_014491</td><td>2AS5_F</td><td>NM_012340</td><td>2AS5_N</td></tr>
<tr><td>123</td><td>NM_014009</td><td>3QRF_F</td><td>NM_012340</td><td>3QRF_N</td></tr>
<tr><td>124</td><td>NM_005479</td><td>5OY4_X</td><td>GSK3B</td><td>5OY4_A</td></tr>
<tr><td>125</td><td>FSHR</td><td>1XWD_C</td><td>NM_000510</td><td>1XWD_B</td></tr>
<tr><td>126</td><td>NM_013365</td><td>1X79_A</td><td>NM_004703</td><td>1X79_B</td></tr>
<tr><td>127</td><td>NM_000515</td><td>1BP3_A</td><td>NM_000949</td><td>1BP3_B</td></tr>
<tr><td>128</td><td>NM_000515</td><td>1BP3_A</td><td>NM_001204317</td><td>1BP3_B</td></tr>
<tr><td>129</td><td>NM_000516</td><td>6NIY_A</td><td>NM_002074</td><td>6NIY_B</td></tr>
<tr><td>130</td><td>NM_001243774</td><td>6NIY_G</td><td>NM_001742</td><td>6NIY_R</td></tr>
<tr><td>131</td><td>NM_001243773</td><td>6NIY_G</td><td>NM_001742</td><td>6NIY_R</td></tr>
<tr><td>132</td><td>NM_053064</td><td>6NIY_G</td><td>NM_001742</td><td>6NIY_R</td></tr>
<tr><td>133</td><td>NM_212492</td><td>4D10_A</td><td>NM_003653</td><td>4D10_C</td></tr>
<tr><td>134</td><td>NM_004489</td><td>2L5G_A</td><td>NM_006312</td><td>2L5G_B</td></tr>
<tr><td>135</td><td>NM_007327</td><td>5H8F_B</td><td>NM_000833</td><td>5H8F_A</td></tr>
<tr><td>136</td><td>NM_007327</td><td>5H8F_B</td><td>NM_001134407</td><td>5H8F_A</td></tr>
<tr><td>137</td><td>NM_000848</td><td>3GTU_A</td><td>NM_000849</td><td>3GTU_B</td></tr>
<tr><td>138</td><td>NM_176795</td><td>1WQ1_R</td><td>NM_002890</td><td>1WQ1_G</td></tr>
<tr><td>139</td><td>NM_176795</td><td>4K81_B</td><td>NM_004490</td><td>4K81_A</td></tr>
<tr><td>140</td><td>NM_176795</td><td>2C5L_A</td><td>NM_016341</td><td>2C5L_C</td></tr>
<tr><td>141</td><td>NM_176795</td><td>4URU_R</td><td>NM_005633</td><td>4URU_S</td></tr>
<tr><td>142</td><td>NM_176795</td><td>6AXG_B</td><td>NM_170604</td><td>6AXG_A</td></tr>
<tr><td>143</td><td>NM_176795</td><td>1HE8_B</td><td>NM_001282427</td><td>1HE8_A</td></tr>
<tr><td>144</td><td>NM_002162</td><td>1T0P_B</td><td>NM_002209</td><td>1T0P_A</td></tr>
<tr><td>145</td><td>NM_000612</td><td>2L29_B</td><td>NM_000876</td><td>2L29_A</td></tr>
<tr><td>146</td><td>NM_000600</td><td>5FUC_A</td><td>NM_000565</td><td>5FUC_C</td></tr>
<tr><td>147</td><td>NM_002184</td><td>1P9M_A</td><td>NM_000600</td><td>1P9M_B</td></tr>
<tr><td>148</td><td>NM_001567</td><td>2KSO_B</td><td>NM_004431</td><td>2KSO_A</td></tr>
<tr><td>149</td><td>NM_014652</td><td>3ZJY_B</td><td>NM_001412</td><td>3ZJY_C</td></tr>
<tr><td>150</td><td>NM_014652</td><td>3ZJY_B</td><td>NM_006325</td><td>3ZJY_A</td></tr>
<tr><td>151</td><td>NM_001571</td><td>1ZOQ_A</td><td>NM_004380</td><td>1ZOQ_C</td></tr>
<tr><td>152</td><td>NM_002227</td><td>5IXD_A</td><td>NM_170743</td><td>5IXD_B</td></tr>
<tr><td>153</td><td>NM_001322195</td><td>6E2Q_A</td><td>NM_000121</td><td>6E2Q_M</td></tr>
<tr><td>154</td><td>NM_001322194</td><td>6E2Q_A</td><td>NM_000121</td><td>6E2Q_M</td></tr>
<tr><td>155</td><td>NM_004972</td><td>6E2Q_A</td><td>NM_000121</td><td>6E2Q_M</td></tr>
<tr><td>156</td><td>NM_001322196</td><td>6E2Q_A</td><td>NM_000121</td><td>6E2Q_M</td></tr>
<tr><td>157</td><td>NM_002253</td><td>2X1X_R</td><td>NM_005429</td><td>2X1X_E</td></tr>
<tr><td>158</td><td>NM_012289</td><td>5NLB_A</td><td>NM_003590</td><td>5NLB_B</td></tr>
<tr><td>159</td><td>NM_000222</td><td>2E9W_A</td><td>NM_003994</td><td>2E9W_C</td></tr>
<tr><td>160</td><td>NM_170606</td><td>5F6K_C</td><td>NM_004674</td><td>5F6K_A</td></tr>
<tr><td>161</td><td>NM_033360</td><td>6EPP_R</td><td>NM_005633</td><td>6EPP_S</td></tr>
<tr><td>162</td><td>NM_004690</td><td>5BRK_B</td><td>NM_018221</td><td>5BRK_A</td></tr>
<tr><td>163</td><td>NM_014572</td><td>4ZRI_C</td><td>NM_000268</td><td>4ZRI_A</td></tr>
<tr><td>164</td><td>NM_018490</td><td>4KT1_A</td><td>RSPO1</td><td>4KT1_E</td></tr>
<tr><td>165</td><td>NM_002370</td><td>2HYI_A</td><td>NM_007359</td><td>2HYI_D</td></tr>
<tr><td>166</td><td>NM_002370</td><td>2HYI_A</td><td>NM_014740</td><td>2HYI_C</td></tr>
<tr><td>167</td><td>NM_002370</td><td>2HYI_A</td><td>NM_005105</td><td>2HYI_B</td></tr>
<tr><td>168</td><td>NM_002745</td><td>4IZ7_C</td><td>NM_003768</td><td>4IZ7_B</td></tr>
<tr><td>169</td><td>NM_138957</td><td>4IZ7_C</td><td>NM_003768</td><td>4IZ7_B</td></tr>
<tr><td>170</td><td>NM_002755</td><td>6U2G_A</td><td>NM_004333</td><td>6U2G_B</td></tr>
<tr><td>171</td><td>NM_002755</td><td>2Y4I_C</td><td>NM_173598</td><td>2Y4I_B</td></tr>
<tr><td>172</td><td>NM_015112</td><td>2KYL_A</td><td>NM_000314</td><td>2KYL_B</td></tr>
<tr><td>173</td><td>NM_197957</td><td>1NLW_B</td><td>NM_002357</td><td>1NLW_A</td></tr>
<tr><td>174</td><td>NM_021960</td><td>3PK1_A</td><td>NM_138763</td><td>3PK1_B</td></tr>
<tr><td>175</td><td>NM_004526</td><td>5JA4_C</td><td>NM_013432</td><td>5JA4_D</td></tr>
<tr><td>176</td><td>MDM2</td><td>1YCR_A</td><td>NM_000546</td><td>1YCR_B</td></tr>
<tr><td>177</td><td>MDM2</td><td>2VJE_A</td><td>NM_002393</td><td>2VJE_B</td></tr>
<tr><td>178</td><td>MDM2</td><td>2MPS_A</td><td>NM_005427</td><td>2MPS_B</td></tr>
<tr><td>179</td><td>NM_130802</td><td>3U88_A</td><td>NM_005933</td><td>3U88_M</td></tr>
<tr><td>180</td><td>NM_000244</td><td>3U88_A</td><td>NM_005933</td><td>3U88_M</td></tr>
<tr><td>181</td><td>NM_130801</td><td>3U88_A</td><td>NM_005933</td><td>3U88_M</td></tr>
<tr><td>182</td><td>NM_130800</td><td>3U88_A</td><td>NM_005933</td><td>3U88_M</td></tr>
<tr><td>183</td><td>NM_130804</td><td>3U88_A</td><td>NM_005933</td><td>3U88_M</td></tr>
<tr><td>184</td><td>NM_130803</td><td>3U88_A</td><td>NM_005933</td><td>3U88_M</td></tr>
<tr><td>185</td><td>NM_000245</td><td>1SHY_B</td><td>NM_000601</td><td>1SHY_A</td></tr>
<tr><td>186</td><td>NM_000179</td><td>2O8B_B</td><td>NM_000251</td><td>2O8B_A</td></tr>
<tr><td>187</td><td>NM_020998</td><td>4QT8_C</td><td>NM_002447</td><td>4QT8_B</td></tr>
<tr><td>188</td><td>MTOR</td><td>5WBU_B</td><td>NM_022372</td><td>5WBU_D</td></tr>
<tr><td>189</td><td>MTOR</td><td>5H64_A</td><td>NM_020761</td><td>5H64_B</td></tr>
<tr><td>190</td><td>MTOR</td><td>2FAP_B</td><td>NM_000801</td><td>2FAP_A</td></tr>
<tr><td>191</td><td>MTOR</td><td>2FAP_B</td><td>NM_054014</td><td>2FAP_A</td></tr>
<tr><td>192</td><td>MTOR</td><td>5WBH_A</td><td>NM_003161</td><td>5WBH_W</td></tr>
<tr><td>193</td><td>NM_002467</td><td>1NKP_A</td><td>NM_197957</td><td>1NKP_B</td></tr>
<tr><td>194</td><td>NM_002468</td><td>3MOP_A</td><td>NM_001570</td><td>3MOP_K</td></tr>
<tr><td>195</td><td>NM_000260</td><td>5MV9_A</td><td>NM_005709</td><td>5MV9_B</td></tr>
<tr><td>196</td><td>NM_022346</td><td>6IGX_B</td><td>NM_015341</td><td>6IGX_A</td></tr>
<tr><td>197</td><td>NM_002486</td><td>1H2T_C</td><td>NM_007362</td><td>1H2T_Z</td></tr>
<tr><td>198</td><td>NM_000265</td><td>1WLP_B</td><td>NM_000101</td><td>1WLP_A</td></tr>
<tr><td>199</td><td>NM_006156</td><td>1XT9_B</td><td>NM_145204</td><td>1XT9_A</td></tr>
<tr><td>200</td><td>NM_005450</td><td>1M4U_A</td><td>NM_001719</td><td>1M4U_L</td></tr>
<tr><td>201</td><td>NM_017617</td><td>2F8X_K</td><td>NM_014757</td><td>2F8X_M</td></tr>
<tr><td>202</td><td>NM_017617</td><td>2F8X_K</td><td>NM_015874</td><td>2F8X_C</td></tr>
<tr><td>203</td><td>NM_017617</td><td>5L0R_B</td><td>NM_152305</td><td>5L0R_A</td></tr>
<tr><td>204</td><td>NM_017617</td><td>6PY8_F</td><td>NM_001004354</td><td>6PY8_B</td></tr>
<tr><td>205</td><td>NM_002527</td><td>1BND_B</td><td>NM_170735</td><td>1BND_A</td></tr>
<tr><td>206</td><td>NM_002529</td><td>1WWW_X</td><td>NM_002506</td><td>1WWW_V</td></tr>
<tr><td>207</td><td>NM_002568</td><td>1JGN_A</td><td>NM_016480</td><td>1JGN_B</td></tr>
<tr><td>208</td><td>NM_016948</td><td>1WMH_B</td><td>NM_002740</td><td>1WMH_A</td></tr>
<tr><td>209</td><td>NM_002585</td><td>1B72_B</td><td>NM_002144</td><td>1B72_A</td></tr>
<tr><td>210</td><td>NM_002609</td><td>3MJG_X</td><td>NM_002608</td><td>3MJG_A</td></tr>
<tr><td>211</td><td>NM_003477</td><td>1ZY8_K</td><td>NM_000108</td><td>1ZY8_A</td></tr>
<tr><td>212</td><td>NM_032758</td><td>5IFE_D</td><td>NM_012433</td><td>5IFE_C</td></tr>
<tr><td>213</td><td>NM_006218</td><td>6NCT_A</td><td>NM_181504</td><td>6NCT_B</td></tr>
<tr><td>214</td><td>NM_005026</td><td>6PYU_A</td><td>NM_181504</td><td>6PYU_B</td></tr>
<tr><td>215</td><td>NM_002669</td><td>5MQF_D</td><td>NM_014691</td><td>5MQF_U</td></tr>
<tr><td>216</td><td>NM_006231</td><td>5VBN_B</td><td>NM_002692</td><td>5VBN_A</td></tr>
<tr><td>217</td><td>NM_005037</td><td>2PRG_A</td><td>NM_147223</td><td>2PRG_C</td></tr>
<tr><td>218</td><td>NM_014225</td><td>2IE3_A</td><td>NM_002715</td><td>2IE3_C</td></tr>
<tr><td>219</td><td>PRPF8</td><td>5MQF_A</td><td>NM_004247</td><td>5MQF_B</td></tr>
<tr><td>220</td><td>NM_000264</td><td>6DMY_A</td><td>NM_000193</td><td>6DMY_B</td></tr>
<tr><td>221</td><td>NM_005051</td><td>4R3Z_C</td><td>NM_002887</td><td>4R3Z_B</td></tr>
<tr><td>222</td><td>NM_018890</td><td>2FJU_A</td><td>NM_004573</td><td>2FJU_B</td></tr>
<tr><td>223</td><td>NM_018890</td><td>2FJU_A</td><td>NM_001284299</td><td>2FJU_B</td></tr>
<tr><td>224</td><td>NM_018890</td><td>4YON_B</td><td>NM_020820</td><td>4YON_A</td></tr>
<tr><td>225</td><td>NM_018890</td><td>3B13_B</td><td>NM_004946</td><td>3B13_A</td></tr>
<tr><td>226</td><td>NM_018890</td><td>3SUA_A</td><td>NM_002673</td><td>3SUA_D</td></tr>
<tr><td>227</td><td>NM_018890</td><td>1E96_A</td><td>NM_000433</td><td>1E96_B</td></tr>
<tr><td>228</td><td>NM_020165</td><td>2YBF_B</td><td>NM_003337</td><td>2YBF_A</td></tr>
<tr><td>229</td><td>NM_020165</td><td>2MRE_B</td><td>NM_021009</td><td>2MRE_A</td></tr>
<tr><td>230</td><td>NM_006265</td><td>4PK7_B</td><td>NM_006603</td><td>4PK7_A</td></tr>
<tr><td>231</td><td>NM_002880</td><td>1C1Y_B</td><td>NM_002884</td><td>1C1Y_A</td></tr>
<tr><td>232</td><td>NM_002880</td><td>4G0N_B</td><td>NM_176795</td><td>4G0N_A</td></tr>
<tr><td>233</td><td>NM_018047</td><td>5MQF_P</td><td>NM_016652</td><td>5MQF_O</td></tr>
<tr><td>234</td><td>NM_000321</td><td>2AZE_C</td><td>NM_007111</td><td>2AZE_A</td></tr>
<tr><td>235</td><td>NM_000321</td><td>2AZE_C</td><td>NM_005225</td><td>2AZE_B</td></tr>
<tr><td>236</td><td>NM_021975</td><td>2LSP_A</td><td>NM_058243</td><td>2LSP_B</td></tr>
<tr><td>237</td><td>NM_020975</td><td>6GL7_F</td><td>NM_001495</td><td>6GL7_D</td></tr>
<tr><td>238</td><td>NM_020975</td><td>6GL7_F</td><td>NM_004558</td><td>6GL7_B</td></tr>
<tr><td>239</td><td>NM_005614</td><td>3T5G_A</td><td>NM_002601</td><td>3T5G_B</td></tr>
<tr><td>240</td><td>NM_001664</td><td>1CXZ_A</td><td>NM_002741</td><td>1CXZ_B</td></tr>
<tr><td>241</td><td>NM_001664</td><td>1CXZ_A</td><td>NM_213560</td><td>1CXZ_B</td></tr>
<tr><td>242</td><td>NM_001664</td><td>6BC0_F</td><td>NM_001177693</td><td>6BC0_A</td></tr>
<tr><td>243</td><td>NM_001664</td><td>4XH9_B</td><td>NM_005863</td><td>4XH9_A</td></tr>
<tr><td>244</td><td>NM_001664</td><td>1OW3_B</td><td>NM_004308</td><td>1OW3_A</td></tr>
<tr><td>245</td><td>NM_001664</td><td>1X86_B</td><td>NM_015313</td><td>1X86_A</td></tr>
<tr><td>246</td><td>NM_017763</td><td>4KNG_E</td><td>NM_003667</td><td>4KNG_A</td></tr>
<tr><td>247</td><td>NM_000328</td><td>4JHP_C</td><td>NM_002601</td><td>4JHP_B</td></tr>
<tr><td>248</td><td>RSPO1</td><td>4KNG_M</td><td>NM_017763</td><td>4KNG_E</td></tr>
<tr><td>249</td><td>NM_001001890</td><td>1H9D_A</td><td>NM_001755</td><td>1H9D_B</td></tr>
<tr><td>250</td><td>NM_002957</td><td>3DZU_A</td><td>NM_005037</td><td>3DZU_D</td></tr>
<tr><td>251</td><td>NM_002957</td><td>4ZO1_B</td><td>NM_000461</td><td>4ZO1_X</td></tr>
<tr><td>252</td><td>NM_021627</td><td>1TGZ_A</td><td>NM_003352</td><td>1TGZ_B</td></tr>
<tr><td>253</td><td>NM_012433</td><td>5IFE_C</td><td>NM_012426</td><td>5IFE_A</td></tr>
<tr><td>254</td><td>NM_012433</td><td>5IFE_C</td><td>NM_031287</td><td>5IFE_B</td></tr>
<tr><td>255</td><td>NM_005414</td><td>5C4V_B</td><td>NM_005359</td><td>5C4V_A</td></tr>
<tr><td>256</td><td>NM_170679</td><td>1P22_B</td><td>NM_033637</td><td>1P22_A</td></tr>
<tr><td>257</td><td>NM_170679</td><td>5JH5_B</td><td>NM_032673</td><td>5JH5_C</td></tr>
<tr><td>258</td><td>NM_005359</td><td>1MR1_A</td><td>NM_003036</td><td>1MR1_C</td></tr>
<tr><td>259</td><td>NM_005359</td><td>1U7V_B</td><td>NM_005901</td><td>1U7V_A</td></tr>
<tr><td>260</td><td>NM_005359</td><td>1U7F_B</td><td>NM_005902</td><td>1U7F_A</td></tr>
<tr><td>261</td><td>NM_003073</td><td>5GJK_B</td><td>NM_003074</td><td>5GJK_A</td></tr>
<tr><td>262</td><td>NM_005631</td><td>6OT0_R</td><td>NM_002069</td><td>6OT0_A</td></tr>
<tr><td>263</td><td>NM_004814</td><td>5MQF_F</td><td>NM_015891</td><td>5MQF_E</td></tr>
<tr><td>264</td><td>NM_003090</td><td>5MQF_W</td><td>NM_016059</td><td>5MQF_V</td></tr>
<tr><td>265</td><td>NM_005701</td><td>3NBY_B</td><td>NM_006325</td><td>3NBY_C</td></tr>
<tr><td>266</td><td>NM_005701</td><td>5DIS_C</td><td>NM_003400</td><td>5DIS_A</td></tr>
<tr><td>267</td><td>NM_012245</td><td>5MQF_C</td><td>NM_002669</td><td>5MQF_D</td></tr>
<tr><td>268</td><td>NM_003106</td><td>1O4X_B</td><td>NM_002697</td><td>1O4X_A</td></tr>
<tr><td>269</td><td>NM_003710</td><td>1YC0_I</td><td>NM_001528</td><td>1YC0_A</td></tr>
<tr><td>270</td><td>NM_003563</td><td>4EOZ_A</td><td>NM_003590</td><td>4EOZ_B</td></tr>
<tr><td>271</td><td>NM_001318537</td><td>5HKY_B</td><td>NM_005188</td><td>5HKY_A</td></tr>
<tr><td>272</td><td>NM_005842</td><td>5HKY_B</td><td>NM_005188</td><td>5HKY_A</td></tr>
<tr><td>273</td><td>NM_001318538</td><td>5HKY_B</td><td>NM_005188</td><td>5HKY_A</td></tr>
<tr><td>274</td><td>NM_001318536</td><td>5HKY_B</td><td>NM_005188</td><td>5HKY_A</td></tr>
<tr><td>275</td><td>NM_003126</td><td>3LBX_A</td><td>NM_001024858</td><td>3LBX_B</td></tr>
<tr><td>276</td><td>NM_003126</td><td>3LBX_A</td><td>NM_001355436</td><td>3LBX_B</td></tr>
<tr><td>277</td><td>NM_016333</td><td>5MQF_S</td><td>NM_020943</td><td>5MQF_T</td></tr>
<tr><td>278</td><td>NM_000455</td><td>2WTK_C</td><td>NM_001003787</td><td>2WTK_B</td></tr>
<tr><td>279</td><td>NM_000455</td><td>2WTK_C</td><td>NM_016289</td><td>2WTK_A</td></tr>
<tr><td>280</td><td>NM_004606</td><td>3AAD_A</td><td>NM_014034</td><td>3AAD_B</td></tr>
<tr><td>281</td><td>NM_016495</td><td>4Z6Y_B</td><td>NM_000368</td><td>4Z6Y_C</td></tr>
<tr><td>282</td><td>NM_003201</td><td>6ERP_C</td><td>NM_005035</td><td>6ERP_A</td></tr>
<tr><td>283</td><td>NM_003242</td><td>3KFD_E</td><td>TGFB1</td><td>3KFD_A</td></tr>
<tr><td>284</td><td>NM_003242</td><td>5TX4_A</td><td>NM_003238</td><td>5TX4_B</td></tr>
<tr><td>285</td><td>NM_003242</td><td>2PJY_B</td><td>NM_003239</td><td>2PJY_A</td></tr>
<tr><td>286</td><td>NM_003254</td><td>3V96_A</td><td>NM_002425</td><td>3V96_B</td></tr>
<tr><td>287</td><td>NM_005781</td><td>1CF4_B</td><td>NM_001791</td><td>1CF4_A</td></tr>
<tr><td>288</td><td>NM_002270</td><td>1QBK_B</td><td>NM_006325</td><td>1QBK_C</td></tr>
<tr><td>289</td><td>NM_002270</td><td>5YVG_A</td><td>NM_004960</td><td>5YVG_X</td></tr>
<tr><td>290</td><td>NM_000546</td><td>2MEJ_B</td><td>NM_138578</td><td>2MEJ_A</td></tr>
<tr><td>291</td><td>NM_000546</td><td>2LY4_B</td><td>NM_002128</td><td>2LY4_A</td></tr>
<tr><td>292</td><td>NM_000546</td><td>2K8F_B</td><td>NM_001429</td><td>2K8F_A</td></tr>
<tr><td>293</td><td>NM_000546</td><td>2RUK_A</td><td>NM_005316</td><td>2RUK_B</td></tr>
<tr><td>294</td><td>NM_000546</td><td>2B3G_B</td><td>NM_002945</td><td>2B3G_A</td></tr>
<tr><td>295</td><td>NM_000546</td><td>1YCS_A</td><td>NM_001031685</td><td>1YCS_B</td></tr>
<tr><td>296</td><td>NM_000546</td><td>1YCS_A</td><td>NM_005426</td><td>1YCS_B</td></tr>
<tr><td>297</td><td>NM_005657</td><td>6CO2_C</td><td>NM_152395</td><td>6CO2_A</td></tr>
<tr><td>298</td><td>NM_021138</td><td>1F3V_B</td><td>NM_001323552</td><td>1F3V_A</td></tr>
<tr><td>299</td><td>NM_021138</td><td>1F3V_B</td><td>NM_003789</td><td>1F3V_A</td></tr>
<tr><td>300</td><td>NM_016399</td><td>6I3V_A</td><td>NM_013237</td><td>6I3V_F</td></tr>
<tr><td>301</td><td>NM_005499</td><td>1Y8Q_B</td><td>NM_005500</td><td>1Y8Q_A</td></tr>
<tr><td>302</td><td>NM_021009</td><td>1XD3_B</td><td>NM_006002</td><td>1XD3_A</td></tr>
<tr><td>303</td><td>UBE2L3</td><td>1FBV_C</td><td>NM_005188</td><td>1FBV_A</td></tr>
<tr><td>304</td><td>NM_003969</td><td>1Y8X_A</td><td>NM_198195</td><td>1Y8X_B</td></tr>
<tr><td>305</td><td>NM_000551</td><td>4WQO_A</td><td>NM_003591</td><td>4WQO_D</td></tr>
<tr><td>306</td><td>NM_000551</td><td>6HR2_B</td><td>NM_003072</td><td>6HR2_A</td></tr>
<tr><td>307</td><td>NM_000551</td><td>1LM8_V</td><td>NM_005648</td><td>1LM8_C</td></tr>
<tr><td>308</td><td>NM_000551</td><td>1LM8_V</td><td>NM_007108</td><td>1LM8_B</td></tr>
<tr><td>309</td><td>NM_000552</td><td>1M10_A</td><td>NM_000173</td><td>1M10_B</td></tr>
<tr><td>310</td><td>NM_018383</td><td>6FBS_B</td><td>NM_006693</td><td>6FBS_C</td></tr>
<tr><td>311</td><td>NM_018383</td><td>6FBS_B</td><td>NM_013291</td><td>6FBS_A</td></tr>
<tr><td>312</td><td>NM_020196</td><td>5MQF_M</td><td>NM_015484</td><td>5MQF_N</td></tr>
<tr><td>313</td><td>NM_003400</td><td>5DIS_A</td><td>NM_006325</td><td>5DIS_B</td></tr>
<tr><td>314</td><td>NM_003400</td><td>2L1L_B</td><td>NM_006823</td><td>2L1L_A</td></tr>
<tr><td>315</td><td>NM_002771</td><td>6BX8_A</td><td>NM_006287</td><td>6BX8_E</td></tr>
<tr><td>316</td><td>NM_002771</td><td>6HAR_A</td><td>NM_000484</td><td>6HAR_E</td></tr>
<tr><td>317</td><td>NM_002771</td><td>5JBT_A</td><td>NM_001642</td><td>5JBT_X</td></tr>
<tr><td>318</td><td>NM_002771</td><td>4U30_A</td><td>NM_001633</td><td>4U30_X</td></tr>
<tr><td>319</td><td>NM_002771</td><td>4U32_A</td><td>SPINT2</td><td>4U32_X</td></tr>
<tr><td>320</td><td>LILRB1</td><td>4NO0_D</td><td>NM_002339</td><td>4NO0_C</td></tr>
<tr><td>321</td><td>LILRB1</td><td>4NO0_D</td><td>NM_004048</td><td>4NO0_B</td></tr>
<tr><td>322</td><td>LILRB2</td><td>2DYP_D</td><td>NM_004048</td><td>2DYP_B</td></tr>
<tr><td>323</td><td>NM_002116</td><td>6MPP_A</td><td>NM_004048</td><td>6MPP_C</td></tr>
<tr><td>324</td><td>NM_002568</td><td>3PTH_A</td><td>NM_015155</td><td>3PTH_B</td></tr>
<tr><td>325</td><td>NM_001080467</td><td>4LWZ_B</td><td>NM_004663</td><td>4LWZ_A</td></tr>
<tr><td>326</td><td>NM_003211</td><td>2D07_A</td><td>NM_006937</td><td>2D07_B</td></tr>
<tr><td>327</td><td>NM_003211</td><td>1WYW_A</td><td>NM_003352</td><td>1WYW_B</td></tr>
<tr><td>328</td><td>NM_013289</td><td>5B38_G</td><td>NM_004048</td><td>5B38_B</td></tr>
<tr><td>329</td><td>NM_003348</td><td>6S53_E</td><td>NM_021009</td><td>6S53_F</td></tr>
<tr><td>330</td><td>NM_002124</td><td>6HBY_B</td><td>NM_019111</td><td>6HBY_A</td></tr>
<tr><td>331</td><td>NM_002417</td><td>2AFF_A</td><td>NM_032390</td><td>2AFF_B</td></tr>
<tr><td>332</td><td>NM_002417</td><td>5J28_C</td><td>NM_002710</td><td>5J28_A</td></tr>
<tr><td>333</td><td>NM_133378</td><td>1YA5_A</td><td>NM_003673</td><td>1YA5_T</td></tr>
<tr><td>334</td><td>NM_133378</td><td>3KNB_A</td><td>NM_015311</td><td>3KNB_B</td></tr>
<tr><td>335</td><td>NM_002956</td><td>3E2U_E</td><td>NM_004082</td><td>3E2U_A</td></tr>
<tr><td>336</td><td>NM_183047</td><td>5Y1Z_C</td><td>NM_080881</td><td>5Y1Z_A</td></tr>
<tr><td>337</td><td>NM_002117</td><td>4NT6_A</td><td>NM_004048</td><td>4NT6_B</td></tr>
<tr><td>338</td><td>NM_004628</td><td>2RVB_A</td><td>NM_005316</td><td>2RVB_B</td></tr>
<tr><td>339</td><td>NM_016320</td><td>3MMY_B</td><td>NM_003610</td><td>3MMY_A</td></tr>
<tr><td>340</td><td>NM_016320</td><td>3MMY_B</td><td>NM_001015885</td><td>3MMY_A</td></tr>
<tr><td>341</td><td>NM_182961</td><td>4DXR_B</td><td>NM_001199579</td><td>4DXR_A</td></tr>
<tr><td>342</td><td>NM_000719</td><td>6DAD_C</td><td>NM_006888</td><td>6DAD_A</td></tr>
<tr><td>343</td><td>NM_004095</td><td>5NVN_B</td><td>NM_004846</td><td>5NVN_A</td></tr>
<tr><td>344</td><td>NM_002332</td><td>2FYL_B</td><td>NM_002337</td><td>2FYL_A</td></tr>
<tr><td>345</td><td>NM_020066</td><td>2YLE_B</td><td>NM_020148</td><td>2YLE_A</td></tr>
<tr><td>346</td><td>NM_006325</td><td>5CLL_A</td><td>NM_006267</td><td>5CLL_B</td></tr>
<tr><td>347</td><td>NM_006267</td><td>1Z5S_D</td><td>NM_002883</td><td>1Z5S_C</td></tr>
<tr><td>348</td><td>NM_138484</td><td>3Q6S_E</td><td>NM_006807</td><td>3Q6S_A</td></tr>
<tr><td>349</td><td>NM_003139</td><td>5L3Q_B</td><td>NM_003136</td><td>5L3Q_A</td></tr>
<tr><td>350</td><td>NM_005154</td><td>3N3K_A</td><td>NM_021009</td><td>3N3K_B</td></tr>
<tr><td>351</td><td>NM_005154</td><td>2GWF_A</td><td>NM_005785</td><td>2GWF_B</td></tr>
<tr><td>352</td><td>NM_000445</td><td>3F7P_A</td><td>ITGB4</td><td>3F7P_C</td></tr>
<tr><td>353</td><td>NM_002156</td><td>4PJ1_A</td><td>NM_002157</td><td>4PJ1_O</td></tr>
<tr><td>354</td><td>NM_031407</td><td>6FYH_A</td><td>NM_018955</td><td>6FYH_B</td></tr>
<tr><td>355</td><td>NM_031407</td><td>5C6H_B</td><td>NM_021960</td><td>5C6H_A</td></tr>
<tr><td>356</td><td>NM_015908</td><td>5OO6_C</td><td>NM_007362</td><td>5OO6_B</td></tr>
<tr><td>357</td><td>NM_006662</td><td>6IGM_H</td><td>NM_022496</td><td>6IGM_G</td></tr>
<tr><td>358</td><td>NM_001103146</td><td>5NVL_B</td><td>NM_004846</td><td>5NVL_A</td></tr>
<tr><td>359</td><td>NM_005652</td><td>3K6G_D</td><td>NM_018975</td><td>3K6G_A</td></tr>
<tr><td>360</td><td>NM_002912</td><td>3VU7_Z</td><td>NM_006341</td><td>3VU7_C</td></tr>
<tr><td>361</td><td>NM_002912</td><td>4GK0_C</td><td>NM_016316</td><td>4GK0_E</td></tr>
<tr><td>362</td><td>NM_001379</td><td>4YOC_A</td><td>NM_003470</td><td>4YOC_C</td></tr>
<tr><td>363</td><td>NM_001379</td><td>5WVO_C</td><td>NM_001177413</td><td>5WVO_A</td></tr>
<tr><td>364</td><td>NM_001379</td><td>5WVO_C</td><td>NM_001135592</td><td>5WVO_A</td></tr>
<tr><td>365</td><td>NM_001379</td><td>5WVO_C</td><td>NM_002954</td><td>5WVO_A</td></tr>
<tr><td>366</td><td>NM_001379</td><td>5YDR_B</td><td>NM_018955</td><td>5YDR_A</td></tr>
<tr><td>367</td><td>NM_152903</td><td>4XC2_E</td><td>NM_007278</td><td>4XC2_A</td></tr>
<tr><td>368</td><td>NM_080703</td><td>4X86_B</td><td>NM_014235</td><td>4X86_A</td></tr>
<tr><td>369</td><td>NM_080703</td><td>6AU8_C</td><td>NM_015949</td><td>6AU8_A</td></tr>
<tr><td>370</td><td>NM_080703</td><td>2N9P_C</td><td>RNF126</td><td>2N9P_A</td></tr>
<tr><td>371</td><td>NM_005085</td><td>3FMP_A</td><td>NM_007242</td><td>3FMP_B</td></tr>
<tr><td>372</td><td>NM_000206</td><td>6OEL_C</td><td>NM_000418</td><td>6OEL_B</td></tr>
<tr><td>373</td><td>NM_000206</td><td>6OEL_C</td><td>NM_001257406</td><td>6OEL_B</td></tr>
<tr><td>374</td><td>NM_000206</td><td>5M5E_C</td><td>NM_000586</td><td>5M5E_D</td></tr>
<tr><td>375</td><td>NM_000206</td><td>5M5E_C</td><td>NM_001346223</td><td>5M5E_B</td></tr>
<tr><td>376</td><td>NM_000206</td><td>5M5E_C</td><td>NM_001346222</td><td>5M5E_B</td></tr>
<tr><td>377</td><td>NM_000206</td><td>5M5E_C</td><td>NM_000878</td><td>5M5E_B</td></tr>
<tr><td>378</td><td>NM_003906</td><td>4DHX_A</td><td>NM_020189</td><td>4DHX_B</td></tr>
<tr><td>379</td><td>NM_001164</td><td>3DXC_A</td><td>NM_000484</td><td>3DXC_B</td></tr>
<tr><td>380</td><td>NM_003758</td><td>2KRB_B</td><td>NM_001037283</td><td>2KRB_A</td></tr>
<tr><td>381</td><td>NM_003758</td><td>2KRB_B</td><td>NM_003751</td><td>2KRB_A</td></tr>
<tr><td>382</td><td>NM_003758</td><td>2KRB_B</td><td>NM_001362791</td><td>2KRB_A</td></tr>
<tr><td>383</td><td>NM_003160</td><td>6GR8_A</td><td>NM_020238</td><td>6GR8_B</td></tr>
<tr><td>384</td><td>NM_178172</td><td>6OAU_C</td><td>NM_000237</td><td>6OAU_A</td></tr>
<tr><td>385</td><td>NM_002123</td><td>6DFX_B</td><td>NM_002122</td><td>6DFX_A</td></tr>
<tr><td>386</td><td>NM_207410</td><td>6Q2J_C</td><td>NM_020975</td><td>6Q2J_E</td></tr>
<tr><td>387</td><td>NM_013444</td><td>6MUN_B</td><td>NM_002810</td><td>6MUN_A</td></tr>
<tr><td>388</td><td>NM_006082</td><td>6I2I_A</td><td>NM_178014</td><td>6I2I_B</td></tr>
<tr><td>389</td><td>NM_004999</td><td>6J56_A</td><td>NM_005488</td><td>6J56_C</td></tr>
<tr><td>390</td><td>NM_002759</td><td>6NPY_B</td><td>NM_004895</td><td>6NPY_A</td></tr>
<tr><td>391</td><td>NM_014985</td><td>6CSU_B</td><td>NM_025180</td><td>6CSU_C</td></tr>
<tr><td>392</td><td>NM_004397</td><td>6S8S_A</td><td>NM_025083</td><td>6S8S_B</td></tr>
<tr><td>393</td><td>NM_005862</td><td>6RRK_B</td><td>NM_006265</td><td>6RRK_C</td></tr>
<tr><td>394</td><td>NM_021960</td><td>6QFI_A</td><td>NM_001204108</td><td>6QFI_B</td></tr>
<tr><td>395</td><td>NM_021960</td><td>6QFI_A</td><td>NM_138621</td><td>6QFI_B</td></tr>
<tr><td>396</td><td>NM_021960</td><td>6QFI_A</td><td>NM_138622</td><td>6QFI_B</td></tr>
<tr><td>397</td><td>NM_021960</td><td>6QFI_A</td><td>NM_138627</td><td>6QFI_B</td></tr>
<tr><td>398</td><td>NM_003907</td><td>6O81_A</td><td>NM_014239</td><td>6O81_C</td></tr>
<tr><td>399</td><td>NM_006663</td><td>6DCX_C</td><td>NM_002708</td><td>6DCX_A</td></tr>
<tr><td>400</td><td>NM_018225</td><td>6Q8I_A</td><td>NM_006083</td><td>6Q8I_C</td></tr>
<tr><td>401</td><td>NM_138357</td><td>6O5B_J</td><td>NM_033318</td><td>6O5B_F</td></tr>
<tr><td>402</td><td>MTREX</td><td>6IEH_B</td><td>NM_017970</td><td>6IEH_A</td></tr>
<tr><td>403</td><td>NM_000602</td><td>5ZLZ_I</td><td>NM_000930</td><td>5ZLZ_E</td></tr>
<tr><td>404</td><td>NM_033407</td><td>6AJ4_A</td><td>NM_001791</td><td>6AJ4_B</td></tr>
<tr><td>405</td><td>NM_001031685</td><td>6GHM_C</td><td>NM_002708</td><td>6GHM_A</td></tr>
<tr><td>406</td><td>NM_005426</td><td>6GHM_C</td><td>NM_002708</td><td>6GHM_A</td></tr>
<tr><td>407</td><td>NM_018489</td><td>6AGO_A</td><td>NM_006791</td><td>6AGO_C</td></tr>
<tr><td>408</td><td>NM_004844</td><td>6IXV_A</td><td>NM_004663</td><td>6IXV_E</td></tr>
<tr><td>409</td><td>NM_002265</td><td>6N89_A</td><td>NM_005318</td><td>6N89_B</td></tr>
<tr><td>410</td><td>NM_005125</td><td>6FP6_B</td><td>NM_000454</td><td>6FP6_A</td></tr>
<tr><td>411</td><td>NM_007327</td><td>6IRA_A</td><td>NM_000833</td><td>6IRA_D</td></tr>
<tr><td>412</td><td>NM_007327</td><td>6IRA_A</td><td>NM_001134407</td><td>6IRA_D</td></tr>
<tr><td>413</td><td>NM_000537</td><td>6I3F_B</td><td>NM_000029</td><td>6I3F_A</td></tr>
<tr><td>414</td><td>NM_031966</td><td>6GU4_B</td><td>NM_001786</td><td>6GU4_A</td></tr>
<tr><td>415</td><td>NM_005901</td><td>5ZOJ_A</td><td>LEMD3</td><td>5ZOJ_D</td></tr>
<tr><td>416</td><td>NM_006861</td><td>6EKK_C</td><td>NM_024820</td><td>6EKK_A</td></tr>
<tr><td>417</td><td>NM_005123</td><td>6A5Z_A</td><td>NM_002957</td><td>6A5Z_D</td></tr>
<tr><td>418</td><td>NM_000334</td><td>6AGF_A</td><td>NM_001037</td><td>6AGF_B</td></tr>
<tr><td>419</td><td>NM_178012</td><td>6E7C_B</td><td>NM_006082</td><td>6E7C_A</td></tr>
<tr><td>420</td><td>NM_000183</td><td>6DV2_A</td><td>NM_000182</td><td>6DV2_G</td></tr>
<tr><td>421</td><td>NM_002502</td><td>5ZMC_A</td><td>NM_005238</td><td>5ZMC_B</td></tr>
<tr><td>422</td><td>NM_001261403</td><td>5ZMC_A</td><td>NM_005238</td><td>5ZMC_B</td></tr>
<tr><td>423</td><td>NM_001288724</td><td>5ZMC_A</td><td>NM_005238</td><td>5ZMC_B</td></tr>
<tr><td>424</td><td>NM_001682</td><td>6A69_A</td><td>NM_012428</td><td>6A69_B</td></tr>
<tr><td>425</td><td>NM_001240</td><td>6GZH_B</td><td>NM_001261</td><td>6GZH_A</td></tr>
<tr><td>426</td><td>NM_004990</td><td>5Y6L_A</td><td>NM_004446</td><td>5Y6L_B</td></tr>
<tr><td>427</td><td>NM_000297</td><td>6A70_A</td><td>NM_001009944</td><td>6A70_B</td></tr>
<tr><td>428</td><td>NM_000090</td><td>6FZW_A</td><td>NM_002593</td><td>6FZW_D</td></tr>
<tr><td>429</td><td>NM_014741</td><td>5XV6_A</td><td>NM_021934</td><td>5XV6_B</td></tr>
<tr><td>430</td><td>NM_021723</td><td>5Y31_A</td><td>NM_005097</td><td>5Y31_B</td></tr>
</table>
