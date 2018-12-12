

```python
%matplotlib inline
from matplotlib.pyplot import figure, imshow, axis
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
import numpy as np
from sklearn.metrics import mean_squared_error
from sklearn.metrics import r2_score
from math import sqrt
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split
from math import log10
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import Draw
from IPython.display import display, Markdown, HTML
```

** Require RDKit 2018.09.01 or later **


```python
import rdkit
rdkit.__version__
```




    '2018.09.1'




```python
#sdf_file = "CHEMBL1827733_5HT2A.sdf"
sdf_file = "CHEMBL930273_GSK3.sdf"
#sdf_file = "CHEMBL952131_EGFR.sdf"
```


```python
def sdf_to_desc(sdf_file):
    fps = []
    targets = []
    nfps = []
    mols=[]
    bis = []

    for mol in Chem.SDMolSupplier(sdf_file):
        mols.append(mol)
        bi = {}
        fps.append(AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, bitInfo=bi))  # fingerprint
        targets.append(9.0 - log10(float(mol.GetProp("ACTIVITY"))))  # pIC50
        bis.append(bi)

    for fp in fps:
        nfp = np.zeros((1,))
        DataStructs.ConvertToNumpyArray(fp, nfp)
        nfps.append(nfp)

    return (np.array(nfps), np.array(targets), mols, bis)
```

### Random Forest Regression Analysis


```python
def rf(x_train, x_test, y_train, y_test):
    r = RandomForestRegressor().fit(x_train, y_train)
    y_pred = r.predict(x_test)
    r2 = r2_score(y_test, y_pred)
    rmse = sqrt(mean_squared_error(y_test, y_pred))
    print(' RF: R2: {0:f}, RMSE:{1:f}'.format(r2, rmse))
    return r
```


```python
x, y, mols,bis = sdf_to_desc(sdf_file)
```


```python
x_train, x_test, y_train, y_test, mols_train, mols_test, bis_train, bis_test = train_test_split(x, y, mols, bis, test_size=0.1)
```


```python
r = rf(x_train, x_test, y_train, y_test)
```

     RF: R2: 0.481589, RMSE:0.802244


### Sorting bit positions by feature importances


```python
sorted_pos = sorted(zip(range(r.feature_importances_.shape[0]), r.feature_importances_), key=lambda x: x[1], reverse=True)
```


```python
important_posisions = [t[0] for t in sorted_pos[:20]]
```

### Visualizing important fragments


```python
i = 0
for bi, mol in zip(bis_test, mols_test):
    i += 1
    display(Markdown("## Test molecule #{}<h2>".format(i)))
    display(mol)
    display(Markdown("#### Important fragments"))
    frgs = []
    for pos in important_posisions:
        if bi.get(pos, False):
            print("Bit position: {}".format(pos))
            display(Draw.DrawMorganBit(mol,pos,bi))

```


## Test molecule #1<h2>



![png](Visualizing_important_fragments_files/Visualizing_important_fragments_14_1.png)



#### Important fragments


    Bit position: 1152



![png](Visualizing_important_fragments_files/Visualizing_important_fragments_14_4.png)


    Bit position: 155



![png](Visualizing_important_fragments_files/Visualizing_important_fragments_14_6.png)


    Bit position: 650



![png](Visualizing_important_fragments_files/Visualizing_important_fragments_14_8.png)



## Test molecule #2<h2>



![png](Visualizing_important_fragments_files/Visualizing_important_fragments_14_10.png)



#### Important fragments


    Bit position: 328



![png](Visualizing_important_fragments_files/Visualizing_important_fragments_14_13.png)


    Bit position: 926



![png](Visualizing_important_fragments_files/Visualizing_important_fragments_14_15.png)


    Bit position: 591



![png](Visualizing_important_fragments_files/Visualizing_important_fragments_14_17.png)


    Bit position: 1490



![png](Visualizing_important_fragments_files/Visualizing_important_fragments_14_19.png)


    Bit position: 1816



![png](Visualizing_important_fragments_files/Visualizing_important_fragments_14_21.png)


    Bit position: 1152



![png](Visualizing_important_fragments_files/Visualizing_important_fragments_14_23.png)


    Bit position: 898



![png](Visualizing_important_fragments_files/Visualizing_important_fragments_14_25.png)


    Bit position: 650



![png](Visualizing_important_fragments_files/Visualizing_important_fragments_14_27.png)


    Bit position: 896



![png](Visualizing_important_fragments_files/Visualizing_important_fragments_14_29.png)


    Bit position: 980



![png](Visualizing_important_fragments_files/Visualizing_important_fragments_14_31.png)



## Test molecule #3<h2>



![png](Visualizing_important_fragments_files/Visualizing_important_fragments_14_33.png)



#### Important fragments


    Bit position: 1602



![png](Visualizing_important_fragments_files/Visualizing_important_fragments_14_36.png)


    Bit position: 1152



![png](Visualizing_important_fragments_files/Visualizing_important_fragments_14_38.png)


    Bit position: 155



![png](Visualizing_important_fragments_files/Visualizing_important_fragments_14_40.png)


    Bit position: 650



![png](Visualizing_important_fragments_files/Visualizing_important_fragments_14_42.png)



## Test molecule #4<h2>



![png](Visualizing_important_fragments_files/Visualizing_important_fragments_14_44.png)



#### Important fragments


    Bit position: 328



![png](Visualizing_important_fragments_files/Visualizing_important_fragments_14_47.png)


    Bit position: 1602



![png](Visualizing_important_fragments_files/Visualizing_important_fragments_14_49.png)


    Bit position: 926



![png](Visualizing_important_fragments_files/Visualizing_important_fragments_14_51.png)


    Bit position: 591



![png](Visualizing_important_fragments_files/Visualizing_important_fragments_14_53.png)


    Bit position: 1313



![png](Visualizing_important_fragments_files/Visualizing_important_fragments_14_55.png)


    Bit position: 1490



![png](Visualizing_important_fragments_files/Visualizing_important_fragments_14_57.png)


    Bit position: 786



![png](Visualizing_important_fragments_files/Visualizing_important_fragments_14_59.png)


    Bit position: 1816



![png](Visualizing_important_fragments_files/Visualizing_important_fragments_14_61.png)


    Bit position: 1152



![png](Visualizing_important_fragments_files/Visualizing_important_fragments_14_63.png)


    Bit position: 898



![png](Visualizing_important_fragments_files/Visualizing_important_fragments_14_65.png)


    Bit position: 650



![png](Visualizing_important_fragments_files/Visualizing_important_fragments_14_67.png)


    Bit position: 1121



![png](Visualizing_important_fragments_files/Visualizing_important_fragments_14_69.png)


    Bit position: 896



![png](Visualizing_important_fragments_files/Visualizing_important_fragments_14_71.png)



## Test molecule #5<h2>



![png](Visualizing_important_fragments_files/Visualizing_important_fragments_14_73.png)



#### Important fragments


    Bit position: 328



![png](Visualizing_important_fragments_files/Visualizing_important_fragments_14_76.png)


    Bit position: 591



![png](Visualizing_important_fragments_files/Visualizing_important_fragments_14_78.png)


    Bit position: 1452



![png](Visualizing_important_fragments_files/Visualizing_important_fragments_14_80.png)


    Bit position: 786



![png](Visualizing_important_fragments_files/Visualizing_important_fragments_14_82.png)


    Bit position: 1816



![png](Visualizing_important_fragments_files/Visualizing_important_fragments_14_84.png)


    Bit position: 1152



![png](Visualizing_important_fragments_files/Visualizing_important_fragments_14_86.png)


    Bit position: 1535



![png](Visualizing_important_fragments_files/Visualizing_important_fragments_14_88.png)


    Bit position: 650



![png](Visualizing_important_fragments_files/Visualizing_important_fragments_14_90.png)



## Test molecule #6<h2>



![png](Visualizing_important_fragments_files/Visualizing_important_fragments_14_92.png)



#### Important fragments


    Bit position: 1452



![png](Visualizing_important_fragments_files/Visualizing_important_fragments_14_95.png)


    Bit position: 1535



![png](Visualizing_important_fragments_files/Visualizing_important_fragments_14_97.png)



## Test molecule #7<h2>



![png](Visualizing_important_fragments_files/Visualizing_important_fragments_14_99.png)



#### Important fragments


    Bit position: 1602



![png](Visualizing_important_fragments_files/Visualizing_important_fragments_14_102.png)


    Bit position: 1313



![png](Visualizing_important_fragments_files/Visualizing_important_fragments_14_104.png)


    Bit position: 1152



![png](Visualizing_important_fragments_files/Visualizing_important_fragments_14_106.png)


    Bit position: 155



![png](Visualizing_important_fragments_files/Visualizing_important_fragments_14_108.png)


    Bit position: 650



![png](Visualizing_important_fragments_files/Visualizing_important_fragments_14_110.png)



## Test molecule #8<h2>



![png](Visualizing_important_fragments_files/Visualizing_important_fragments_14_112.png)



#### Important fragments


    Bit position: 926



![png](Visualizing_important_fragments_files/Visualizing_important_fragments_14_115.png)


    Bit position: 650



![png](Visualizing_important_fragments_files/Visualizing_important_fragments_14_117.png)


    Bit position: 896



![png](Visualizing_important_fragments_files/Visualizing_important_fragments_14_119.png)



## Test molecule #9<h2>



![png](Visualizing_important_fragments_files/Visualizing_important_fragments_14_121.png)



#### Important fragments


    Bit position: 926



![png](Visualizing_important_fragments_files/Visualizing_important_fragments_14_124.png)


    Bit position: 650



![png](Visualizing_important_fragments_files/Visualizing_important_fragments_14_126.png)



## Test molecule #10<h2>



![png](Visualizing_important_fragments_files/Visualizing_important_fragments_14_128.png)



#### Important fragments


    Bit position: 328



![png](Visualizing_important_fragments_files/Visualizing_important_fragments_14_131.png)


    Bit position: 926



![png](Visualizing_important_fragments_files/Visualizing_important_fragments_14_133.png)


    Bit position: 591



![png](Visualizing_important_fragments_files/Visualizing_important_fragments_14_135.png)


    Bit position: 1490



![png](Visualizing_important_fragments_files/Visualizing_important_fragments_14_137.png)


    Bit position: 1816



![png](Visualizing_important_fragments_files/Visualizing_important_fragments_14_139.png)


    Bit position: 1152



![png](Visualizing_important_fragments_files/Visualizing_important_fragments_14_141.png)


    Bit position: 898



![png](Visualizing_important_fragments_files/Visualizing_important_fragments_14_143.png)


    Bit position: 650



![png](Visualizing_important_fragments_files/Visualizing_important_fragments_14_145.png)


    Bit position: 980



![png](Visualizing_important_fragments_files/Visualizing_important_fragments_14_147.png)



## Test molecule #11<h2>



![png](Visualizing_important_fragments_files/Visualizing_important_fragments_14_149.png)



#### Important fragments


    Bit position: 926



![png](Visualizing_important_fragments_files/Visualizing_important_fragments_14_152.png)


    Bit position: 366



![png](Visualizing_important_fragments_files/Visualizing_important_fragments_14_154.png)


    Bit position: 650



![png](Visualizing_important_fragments_files/Visualizing_important_fragments_14_156.png)



## Test molecule #12<h2>



![png](Visualizing_important_fragments_files/Visualizing_important_fragments_14_158.png)



#### Important fragments


    Bit position: 1152



![png](Visualizing_important_fragments_files/Visualizing_important_fragments_14_161.png)


    Bit position: 155



![png](Visualizing_important_fragments_files/Visualizing_important_fragments_14_163.png)


    Bit position: 650



![png](Visualizing_important_fragments_files/Visualizing_important_fragments_14_165.png)


    Bit position: 896



![png](Visualizing_important_fragments_files/Visualizing_important_fragments_14_167.png)



## Test molecule #13<h2>



![png](Visualizing_important_fragments_files/Visualizing_important_fragments_14_169.png)



#### Important fragments


    Bit position: 328



![png](Visualizing_important_fragments_files/Visualizing_important_fragments_14_172.png)


    Bit position: 926



![png](Visualizing_important_fragments_files/Visualizing_important_fragments_14_174.png)


    Bit position: 591



![png](Visualizing_important_fragments_files/Visualizing_important_fragments_14_176.png)


    Bit position: 1490



![png](Visualizing_important_fragments_files/Visualizing_important_fragments_14_178.png)


    Bit position: 1816



![png](Visualizing_important_fragments_files/Visualizing_important_fragments_14_180.png)


    Bit position: 1152



![png](Visualizing_important_fragments_files/Visualizing_important_fragments_14_182.png)


    Bit position: 898



![png](Visualizing_important_fragments_files/Visualizing_important_fragments_14_184.png)


    Bit position: 650



![png](Visualizing_important_fragments_files/Visualizing_important_fragments_14_186.png)


    Bit position: 980



![png](Visualizing_important_fragments_files/Visualizing_important_fragments_14_188.png)



## Test molecule #14<h2>



![png](Visualizing_important_fragments_files/Visualizing_important_fragments_14_190.png)



#### Important fragments


    Bit position: 1602



![png](Visualizing_important_fragments_files/Visualizing_important_fragments_14_193.png)


    Bit position: 366



![png](Visualizing_important_fragments_files/Visualizing_important_fragments_14_195.png)


    Bit position: 1152



![png](Visualizing_important_fragments_files/Visualizing_important_fragments_14_197.png)


    Bit position: 155



![png](Visualizing_important_fragments_files/Visualizing_important_fragments_14_199.png)


    Bit position: 650



![png](Visualizing_important_fragments_files/Visualizing_important_fragments_14_201.png)



## Test molecule #15<h2>



![png](Visualizing_important_fragments_files/Visualizing_important_fragments_14_203.png)



#### Important fragments


    Bit position: 328



![png](Visualizing_important_fragments_files/Visualizing_important_fragments_14_206.png)


    Bit position: 926



![png](Visualizing_important_fragments_files/Visualizing_important_fragments_14_208.png)


    Bit position: 591



![png](Visualizing_important_fragments_files/Visualizing_important_fragments_14_210.png)


    Bit position: 1490



![png](Visualizing_important_fragments_files/Visualizing_important_fragments_14_212.png)


    Bit position: 786



![png](Visualizing_important_fragments_files/Visualizing_important_fragments_14_214.png)


    Bit position: 1816



![png](Visualizing_important_fragments_files/Visualizing_important_fragments_14_216.png)


    Bit position: 1152



![png](Visualizing_important_fragments_files/Visualizing_important_fragments_14_218.png)


    Bit position: 650



![png](Visualizing_important_fragments_files/Visualizing_important_fragments_14_220.png)


    Bit position: 1121



![png](Visualizing_important_fragments_files/Visualizing_important_fragments_14_222.png)



## Test molecule #16<h2>



![png](Visualizing_important_fragments_files/Visualizing_important_fragments_14_224.png)



#### Important fragments


    Bit position: 328



![png](Visualizing_important_fragments_files/Visualizing_important_fragments_14_227.png)


    Bit position: 926



![png](Visualizing_important_fragments_files/Visualizing_important_fragments_14_229.png)


    Bit position: 591



![png](Visualizing_important_fragments_files/Visualizing_important_fragments_14_231.png)


    Bit position: 1490



![png](Visualizing_important_fragments_files/Visualizing_important_fragments_14_233.png)


    Bit position: 786



![png](Visualizing_important_fragments_files/Visualizing_important_fragments_14_235.png)


    Bit position: 1816



![png](Visualizing_important_fragments_files/Visualizing_important_fragments_14_237.png)


    Bit position: 1152



![png](Visualizing_important_fragments_files/Visualizing_important_fragments_14_239.png)


    Bit position: 650



![png](Visualizing_important_fragments_files/Visualizing_important_fragments_14_241.png)


    Bit position: 1121



![png](Visualizing_important_fragments_files/Visualizing_important_fragments_14_243.png)

