---
typora-root-url: ./
---

# CancerSRT

​		This repository contains the R and Python script source files involved in the CancerSRT paper, which are used to process and analyze spatial transcriptome data.

## Requirements

R package:

- R 4.2.2
- Seurat >= 4.4.0 and <5.0.x
- SingleR
- CellChat
- copykat
- Other required packages are listed in analyze.R

Python package:

- Python >=3.10
- Pytorch>1.13.1
- scanpy
- numpy
- cell2location
- anndata
- pandas

## Installation

```
git clone https://github.com/qdddddasd/CancerSRT.git  

cd CancerSRT
```

## Scripts Overview

​		The code content is divided into two parts: analysis process and related data.

​		The analysis process code is divided into two parts: R language and Python language. Only the cell2location analysis uses Python language. In the cell2location package, the analysis code is stored in the cell2location.ipynb file, and the results package stores the model generated during the analysis process.

​		The project mainly uses R language. There are three files in the R_analyze package, cellchatv2 analysis and spark analysis are two separate files. The remaining analyses are all saved in the analyze.R file.

​		The related_data folder contains the files needed for analysis. The ensemble.csv file contains the correspondence between the two gene naming methods, which is used to solve the problems that arise in the process of customizing the reference set. The immune.csv file contains the correspondence between cell types and makers, which is used for immune analysis. The signatures.csv file is used for functional status analysis.

## Usage

Introduce how to use the code and precautions during operation

1. Install the required R or Python packages.
2. Select data reading method according to different data formats.
3. Replace the path in the code with the path of your personal data.
4. Run the code of the specified analysis module.

Precautions：

The Seurat version must be SeuratV4, and the subsequent package installation cannot update Seurat. If SeuratV5 is used, the code content needs to be updated, which will affect the analysis process.

## Data

Input: This analysis is based on Seurat objects. The two most required input data for all data analysis processes are expression matrix and position information. Regardless of the data form, these two input data must be used as the basis to build objects for subsequent analysis.

Output: More analysis results. Please see the example for the final result.

![](/result.png)

## Contributing

Contributions are welcome! If you have suggestions or improvements, please open an issue or submit a pull request.

## License

MIT License

Copyright (c) 2024 [Yuying Huo]

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.