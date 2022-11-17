# 【例題2.4】 漸変不等流の一次元解析法

- 標準逐次近似法，直接逐次近似法，Runge-Kutta-Gill（ルンゲ-クッタ-ギル）法による漸変不等流の一次元解析を実施するプログラムです．なお，全てのプログラムは常流のみを対象にしています．

- 計算条件の設定は"inputData.xlsx" で実施します．
  - 計算条件の設定は，”計算条件”シートで実施．
  - 水路形状の設定は，”水路形状”シートで実施．

- 計算は，以下のコマンドで実行できます．

   ~~~ sh
    python 2_4Prog.py
    ~~~

- 計算の実行には，以下のライブラリが必要です．
  - pandas
  - numpy
  - matplotlib
  - scipy
  - openpyxl
  - xlrd
です．
[pipやconda](https://www.python.jp/install/anaconda/pip_and_conda.html)でインストールしてください．

- *.ipynbのファイルは，[jupyter notebook](https://jupyter.org)で実行可能です．
[pipやconda](https://www.python.jp/install/anaconda/pip_and_conda.html)でインストール可能です．また，[VSCode](https://code.visualstudio.com/docs/datascience/jupyter-notebooks)からも使用できます．
