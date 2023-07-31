import pandas as pd

# Test data was randomly generated
MOCK_FILE_1_CONTENT = """@ABC1234567.1 1/1
GCGTGTAATGTTATGATCTTATGCTTGTTTTAGTCCGCTAGGTTCTTTGGTGTACTGCCACTTTTCGATGCCATGCGCATTCTTGGGACTAGGAAGTACGA
+
A*}?_]?+i%(Bh48IfmITF;Wf!@i:q.$$G~H[]lo)d5/Xl;m"I*UDY@uy]_{]2}QFM%2uAd4-K>2{tF,.kVWjcE-f*f?57?)h_flau
@ABC1234567.2 2/1
ACCGCCACGCTTACCGTTTTGGCGCTATGCTTCCATTCTGTTGTCTACGAGGCGATAACAACACGATACGCTCTGTTCTTACTCAGACTTATTCCGAAGCC
+
EQn;SjL`t$a]&GE%S).<9-Hl!#xtwmup/r!LJPP%2b_@SXtaHc'6m>6AmeR#pP~P$;FL1m0XYvI)zKr%2-XlWO\\Y8c.OW!v2./N^\\
@ABC1234567.3 3/1
GGCCCTCGAGATACGCGCGGGAGTACGCTCCCAACTGTTTTATACCCCTGTTTTCATCTTATAGCAACACCCCGGAGTAAGCGCACTCATCCTCTCCTATC
+
h=[U[nDhr,8VmiVsx[F=gZ>UIk5>^4Z/0=N}%vA]U$|{byQAb}:YUH~FJn_92*-Uv$mk*kg,p}SV3>zT/pIK3o-)v[B/%3XcSH&'<
@ABC1234567.4 4/1
GCGTGTAATGTTATGATCTTATGCTTGTTTTAGTCCGCTAGGTTCTTTGGTGTACTGCCACTTTTCGATGCCATGCGCATTCTTGGGACTAGGAAGTACGA
+
'~~3QhIkE}yGC<jzj.Lx1|,oV/{dQMZ/ksH>/iH5uK`3c;K]=a-O/N{3.G+ym=\\_8!-;A`JBA:~-tYt7l]5J3%ZMlf\\},a"NcE}I$
@ABC1234567.5 5/1
ACCGCCACGCATCCTACCTTGTAAGAGGATATCAATGGCGATCGGTGTACAAACAGAGCTGATGCCCACTATTTCACGTAAGTAGTGGGAGGGTCGCGTGC
+
<zpLIB(|@OT;0R`<5n5hbT5<|kSBZo!5kI1)1aL~7qX3\\2MIbj*/'+}pk)kVQdlcK1\\Fe|M68b2`2Y/Zz}"NUhOLzay097Y.$'`]/
"""

MOCK_FILE_2_CONTENT = """@ABC1234567.1 1/2
CTTAAGTCATCGTGAAGGGAGGCATACGTGTTGGACAGAGCCCGAGACGCGGAATTTGCAAGATTTTATGGAAGAGGGGCGGCTCATTGCAAATCGAGCTC
+
mp^]Ro,hU^,R1[\\)2BRU)9,9ksW"UCqAHjE3_W#0*9:`!X.GkWpO{*1^+MLkLRz<gh6c_#db*Olo9j.!2BhNPbjPhVC)l9]hSs#E8
@ABC1234567.2 2/2
AGGCTGCACGGTTTCCGCCCAGGGGATCCGCTTTCGGGTATTTCTCTCACGGTATTTCTCCAAAAAAGAAAATTCTCGTCAGACGCTGGTCAGGATCCAGA
+
goP*&{}8<FDW\\Hk]K2{eAAF{%%:nP}6Woh@"gPqW8`kojH2CciLy0*Vxix7&Zj3V1D"O[A&Uv-v(q,6n-$1dZ6E)#wXiF}!dp")W'
@ABC1234567.3 3/2
TCCAGTACGTACGATTTAGTGAATAGAGCAATCTCAGCGGACTTAAAGAGCGCAAGGCCCATCTTGTATCCATGTTGGTAGGAAGTAATGCACCAGGTGCT
+
52tSeIApSGDrD.gjaVT!'xC([GIVHybb@'OKIc9#&4u%Bp*w;H~JCcm:s,v[vddko@`T'F1,D5X,\\ue.xL6JY=qX>4)Fu?M4+DU7@
@ABC1234567.4 4/2
CACCGCGAATTACAATGCTTGTCTGGGGCAGACGCATTACCAGCTGCNAGAGTCCAGAGTTAAGTTGTGTACCGTTGCCCGTTGGTAGAACTCGCACAAGC
+
{~I&KDf`}L$.=XLjiI)]o/-t5nn$]?mOEB&+`+3-QqMwrG/Y5Jt}X!Gz+b``KzJ!0eJ7#4Y;vnT(>bMfLrwRrlK%gDzh"",0XkU*_
@ABC1234567.5 5/2
GATAGGAGAGGATGAGTGCGCTTACTCCGGGGTGTTGCTATAAGATGAAAACAGGGGTATAAAACAGTTGGGAGCGTACTCCCGCGCGTATCTCGAGGGCC
+
.6qJoOJ-^"e~@U=ZU#gZIYj^Mlp{r)X*$+]"{mUbGl16N(Z,:=4P!09=Z3s5w<tj7V|(9l7,7-E-\\nZ$IIGt*]:)j5R@>/7quD[2O
"""

EXPECTED_RESULT_DF_SINGLE_END = pd.DataFrame(
    {
        "ids": [
            ["@ABC1234567.1 1/1", "@ABC1234567.4 4/1"],
            ["@ABC1234567.2 2/1"],
            ["@ABC1234567.3 3/1"],
            ["@ABC1234567.5 5/1"],
        ],
        "sequence": [
            (
                "GCGTGTAATGTTATGATCTTATGCTTGTTTTAGTCCGCTAGGTTCTTTG"
                "GTGTACTGCCACTTTTCGATGCCATGCGCATTCTTGGGACTAGGAAGTACGA"
            ),
            (
                "ACCGCCACGCTTACCGTTTTGGCGCTATGCTTCCATTCTGTTGTCTACGAG"
                "GCGATAACAACACGATACGCTCTGTTCTTACTCAGACTTATTCCGAAGCC"
            ),
            (
                "GGCCCTCGAGATACGCGCGGGAGTACGCTCCCAACTGTTTTATACCCCTGTT"
                "TTCATCTTATAGCAACACCCCGGAGTAAGCGCACTCATCCTCTCCTATC"
            ),
            (
                "ACCGCCACGCATCCTACCTTGTAAGAGGATATCAATGGCGATCGGTGTACAAAC"
                "AGAGCTGATGCCCACTATTTCACGTAAGTAGTGGGAGGGTCGCGTGC"
            ),
        ],
        "count": [2, 1, 1, 1],
    }
).astype(dtype={"ids": "object", "sequence": "string", "count": "uint32"})

EXPECTED_RESULT_DF_PAIRED_END = pd.DataFrame(
    {
        "ids": [
            ["@ABC1234567.1 1/1", "@ABC1234567.4 4/1"],
            ["@ABC1234567.2 2/1"],
            ["@ABC1234567.3 3/1", "@ABC1234567.5 5/2"],
            ["@ABC1234567.5 5/1"],
            ["@ABC1234567.1 1/2"],
            ["@ABC1234567.2 2/2"],
            ["@ABC1234567.3 3/2"],
            ["@ABC1234567.4 4/2"],
        ],
        "sequence": [
            (
                "GCGTGTAATGTTATGATCTTATGCTTGTTTTAGTCCGCTAGGTTCTTTG"
                "GTGTACTGCCACTTTTCGATGCCATGCGCATTCTTGGGACTAGGAAGTACGA"
            ),
            (
                "ACCGCCACGCTTACCGTTTTGGCGCTATGCTTCCATTCTGTTGTCTACGAG"
                "GCGATAACAACACGATACGCTCTGTTCTTACTCAGACTTATTCCGAAGCC"
            ),
            (
                "GGCCCTCGAGATACGCGCGGGAGTACGCTCCCAACTGTTTTATACCCCTGTT"
                "TTCATCTTATAGCAACACCCCGGAGTAAGCGCACTCATCCTCTCCTATC"
            ),
            (
                "ACCGCCACGCATCCTACCTTGTAAGAGGATATCAATGGCGATCGGTGTACAAAC"
                "AGAGCTGATGCCCACTATTTCACGTAAGTAGTGGGAGGGTCGCGTGC"
            ),
            (
                "GAGCTCGATTTGCAATGAGCCGCCCCTCTTCCATAAAATCTTGCAAATTC"
                "CGCGTCTCGGGCTCTGTCCAACACGTATGCCTCCCTTCACGATGACTTAAG"
            ),
            (
                "TCTGGATCCTGACCAGCGTCTGACGAGAATTTTCTTTTTTGGAGAAATAC"
                "CGTGAGAGAAATACCCGAAAGCGGATCCCCTGGGCGGAAACCGTGCAGCCT"
            ),
            (
                "AGCACCTGGTGCATTACTTCCTACCAACATGGATACAAGATGGGCCTTGCG"
                "CTCTTTAAGTCCGCTGAGATTGCTCTATTCACTAAATCGTACGTACTGGA"
            ),
            (
                "GCTTGTGCGAGTTCTACCAACGGGCAACGGTACACAACTTAACTCTGGACTCT"
                "NGCAGCTGGTAATGCGTCTGCCCCAGACAAGCATTGTAATTCGCGGTG"
            ),
        ],
        "count": [2, 1, 2, 1, 1, 1, 1, 1],
    }
).astype(dtype={"ids": "object", "sequence": "string", "count": "uint32"})

EXPECTED_RESULT_DF_SINGLE_END_CUTOFF = pd.DataFrame(
    {
        "ids": [
            ["@ABC1234567.1 1/1", "@ABC1234567.4 4/1"],
            ["@ABC1234567.2 2/1", "@ABC1234567.5 5/1"],
            ["@ABC1234567.3 3/1"],
        ],
        "sequence": ["GCGTGTAATG", "ACCGCCACGC", "GGCCCTCGAG"],
        "count": [2, 2, 1],
    }
).astype(dtype={"ids": "object", "sequence": "string", "count": "uint32"})

EXPECTED_RESULT_INDEX_SINGLE_END = {
    "ACNVMI": [("@ABC1234567.1 1/1", 0, 0), ("@ABC1234567.4 4/1", 0, 0)],
    "CNVMIL": [("@ABC1234567.1 1/1", 0, 1), ("@ABC1234567.4 4/1", 0, 1)],
    "NVMILC": [("@ABC1234567.1 1/1", 0, 2), ("@ABC1234567.4 4/1", 0, 2)],
    "VMILCL": [("@ABC1234567.1 1/1", 0, 3), ("@ABC1234567.4 4/1", 0, 3)],
    "MILCLF": [("@ABC1234567.1 1/1", 0, 4), ("@ABC1234567.4 4/1", 0, 4)],
    "ILCLF*": [("@ABC1234567.1 1/1", 0, 5), ("@ABC1234567.4 4/1", 0, 5)],
    "LCLF*S": [("@ABC1234567.1 1/1", 0, 6), ("@ABC1234567.4 4/1", 0, 6)],
    "CLF*SA": [("@ABC1234567.1 1/1", 0, 7), ("@ABC1234567.4 4/1", 0, 7)],
    "LF*SAR": [("@ABC1234567.1 1/1", 0, 8), ("@ABC1234567.4 4/1", 0, 8)],
    "F*SARF": [("@ABC1234567.1 1/1", 0, 9), ("@ABC1234567.4 4/1", 0, 9)],
    "*SARFF": [("@ABC1234567.1 1/1", 0, 10), ("@ABC1234567.4 4/1", 0, 10)],
    "SARFFG": [("@ABC1234567.1 1/1", 0, 11), ("@ABC1234567.4 4/1", 0, 11)],
    "ARFFGV": [("@ABC1234567.1 1/1", 0, 12), ("@ABC1234567.4 4/1", 0, 12)],
    "RFFGVL": [("@ABC1234567.1 1/1", 0, 13), ("@ABC1234567.4 4/1", 0, 13)],
    "FFGVLP": [("@ABC1234567.1 1/1", 0, 14), ("@ABC1234567.4 4/1", 0, 14)],
    "FGVLPL": [("@ABC1234567.1 1/1", 0, 15), ("@ABC1234567.4 4/1", 0, 15)],
    "GVLPLF": [("@ABC1234567.1 1/1", 0, 16), ("@ABC1234567.4 4/1", 0, 16)],
    "VLPLFD": [("@ABC1234567.1 1/1", 0, 17), ("@ABC1234567.4 4/1", 0, 17)],
    "LPLFDA": [("@ABC1234567.1 1/1", 0, 18), ("@ABC1234567.4 4/1", 0, 18)],
    "PLFDAM": [("@ABC1234567.1 1/1", 0, 19), ("@ABC1234567.4 4/1", 0, 19)],
    "LFDAMR": [("@ABC1234567.1 1/1", 0, 20), ("@ABC1234567.4 4/1", 0, 20)],
    "FDAMRI": [("@ABC1234567.1 1/1", 0, 21), ("@ABC1234567.4 4/1", 0, 21)],
    "DAMRIL": [("@ABC1234567.1 1/1", 0, 22), ("@ABC1234567.4 4/1", 0, 22)],
    "AMRILG": [("@ABC1234567.1 1/1", 0, 23), ("@ABC1234567.4 4/1", 0, 23)],
    "MRILGT": [("@ABC1234567.1 1/1", 0, 24), ("@ABC1234567.4 4/1", 0, 24)],
    "RILGTR": [("@ABC1234567.1 1/1", 0, 25), ("@ABC1234567.4 4/1", 0, 25)],
    "ILGTRK": [("@ABC1234567.1 1/1", 0, 26), ("@ABC1234567.4 4/1", 0, 26)],
    "LGTRKY": [("@ABC1234567.1 1/1", 0, 27), ("@ABC1234567.4 4/1", 0, 27)],
    "RVML*S": [("@ABC1234567.1 1/1", 1, 0), ("@ABC1234567.4 4/1", 1, 0)],
    "VML*SY": [("@ABC1234567.1 1/1", 1, 1), ("@ABC1234567.4 4/1", 1, 1)],
    "ML*SYA": [("@ABC1234567.1 1/1", 1, 2), ("@ABC1234567.4 4/1", 1, 2)],
    "L*SYAC": [("@ABC1234567.1 1/1", 1, 3), ("@ABC1234567.4 4/1", 1, 3)],
    "*SYACF": [("@ABC1234567.1 1/1", 1, 4), ("@ABC1234567.4 4/1", 1, 4)],
    "SYACFS": [("@ABC1234567.1 1/1", 1, 5), ("@ABC1234567.4 4/1", 1, 5)],
    "YACFSP": [("@ABC1234567.1 1/1", 1, 6), ("@ABC1234567.4 4/1", 1, 6)],
    "ACFSPL": [("@ABC1234567.1 1/1", 1, 7), ("@ABC1234567.4 4/1", 1, 7)],
    "CFSPLG": [("@ABC1234567.1 1/1", 1, 8), ("@ABC1234567.4 4/1", 1, 8)],
    "FSPLGS": [("@ABC1234567.1 1/1", 1, 9), ("@ABC1234567.4 4/1", 1, 9)],
    "SPLGSL": [("@ABC1234567.1 1/1", 1, 10), ("@ABC1234567.4 4/1", 1, 10)],
    "PLGSLV": [("@ABC1234567.1 1/1", 1, 11), ("@ABC1234567.4 4/1", 1, 11)],
    "LGSLVY": [("@ABC1234567.1 1/1", 1, 12), ("@ABC1234567.4 4/1", 1, 12)],
    "GSLVYC": [("@ABC1234567.1 1/1", 1, 13), ("@ABC1234567.4 4/1", 1, 13)],
    "SLVYCH": [("@ABC1234567.1 1/1", 1, 14), ("@ABC1234567.4 4/1", 1, 14)],
    "LVYCHF": [("@ABC1234567.1 1/1", 1, 15), ("@ABC1234567.4 4/1", 1, 15)],
    "VYCHFS": [("@ABC1234567.1 1/1", 1, 16), ("@ABC1234567.4 4/1", 1, 16)],
    "YCHFSM": [("@ABC1234567.1 1/1", 1, 17), ("@ABC1234567.4 4/1", 1, 17)],
    "CHFSMP": [("@ABC1234567.1 1/1", 1, 18), ("@ABC1234567.4 4/1", 1, 18)],
    "HFSMPC": [("@ABC1234567.1 1/1", 1, 19), ("@ABC1234567.4 4/1", 1, 19)],
    "FSMPCA": [("@ABC1234567.1 1/1", 1, 20), ("@ABC1234567.4 4/1", 1, 20)],
    "SMPCAF": [("@ABC1234567.1 1/1", 1, 21), ("@ABC1234567.4 4/1", 1, 21)],
    "MPCAFL": [("@ABC1234567.1 1/1", 1, 22), ("@ABC1234567.4 4/1", 1, 22)],
    "PCAFLG": [("@ABC1234567.1 1/1", 1, 23), ("@ABC1234567.4 4/1", 1, 23)],
    "CAFLGL": [("@ABC1234567.1 1/1", 1, 24), ("@ABC1234567.4 4/1", 1, 24)],
    "AFLGLG": [("@ABC1234567.1 1/1", 1, 25), ("@ABC1234567.4 4/1", 1, 25)],
    "FLGLGS": [("@ABC1234567.1 1/1", 1, 26), ("@ABC1234567.4 4/1", 1, 26)],
    "LGLGST": [("@ABC1234567.1 1/1", 1, 27), ("@ABC1234567.4 4/1", 1, 27)],
    "V*CYDL": [("@ABC1234567.1 1/1", 2, 0), ("@ABC1234567.4 4/1", 2, 0)],
    "*CYDLM": [("@ABC1234567.1 1/1", 2, 1), ("@ABC1234567.4 4/1", 2, 1)],
    "CYDLML": [("@ABC1234567.1 1/1", 2, 2), ("@ABC1234567.4 4/1", 2, 2)],
    "YDLMLV": [("@ABC1234567.1 1/1", 2, 3), ("@ABC1234567.4 4/1", 2, 3)],
    "DLMLVL": [("@ABC1234567.1 1/1", 2, 4), ("@ABC1234567.4 4/1", 2, 4)],
    "LMLVLV": [("@ABC1234567.1 1/1", 2, 5), ("@ABC1234567.4 4/1", 2, 5)],
    "MLVLVR": [("@ABC1234567.1 1/1", 2, 6), ("@ABC1234567.4 4/1", 2, 6)],
    "LVLVR*": [("@ABC1234567.1 1/1", 2, 7), ("@ABC1234567.4 4/1", 2, 7)],
    "VLVR*V": [("@ABC1234567.1 1/1", 2, 8), ("@ABC1234567.4 4/1", 2, 8)],
    "LVR*VL": [("@ABC1234567.1 1/1", 2, 9), ("@ABC1234567.4 4/1", 2, 9)],
    "VR*VLW": [("@ABC1234567.1 1/1", 2, 10), ("@ABC1234567.4 4/1", 2, 10)],
    "R*VLWC": [("@ABC1234567.1 1/1", 2, 11), ("@ABC1234567.4 4/1", 2, 11)],
    "*VLWCT": [("@ABC1234567.1 1/1", 2, 12), ("@ABC1234567.4 4/1", 2, 12)],
    "VLWCTA": [("@ABC1234567.1 1/1", 2, 13), ("@ABC1234567.4 4/1", 2, 13)],
    "LWCTAT": [("@ABC1234567.1 1/1", 2, 14), ("@ABC1234567.4 4/1", 2, 14)],
    "WCTATF": [("@ABC1234567.1 1/1", 2, 15), ("@ABC1234567.4 4/1", 2, 15)],
    "CTATFR": [("@ABC1234567.1 1/1", 2, 16), ("@ABC1234567.4 4/1", 2, 16)],
    "TATFRC": [("@ABC1234567.1 1/1", 2, 17), ("@ABC1234567.4 4/1", 2, 17)],
    "ATFRCH": [("@ABC1234567.1 1/1", 2, 18), ("@ABC1234567.4 4/1", 2, 18)],
    "TFRCHA": [("@ABC1234567.1 1/1", 2, 19), ("@ABC1234567.4 4/1", 2, 19)],
    "FRCHAH": [("@ABC1234567.1 1/1", 2, 20), ("@ABC1234567.4 4/1", 2, 20)],
    "RCHAHS": [("@ABC1234567.1 1/1", 2, 21), ("@ABC1234567.4 4/1", 2, 21)],
    "CHAHSW": [("@ABC1234567.1 1/1", 2, 22), ("@ABC1234567.4 4/1", 2, 22)],
    "HAHSWD": [("@ABC1234567.1 1/1", 2, 23), ("@ABC1234567.4 4/1", 2, 23)],
    "AHSWD*": [("@ABC1234567.1 1/1", 2, 24), ("@ABC1234567.4 4/1", 2, 24)],
    "HSWD*E": [("@ABC1234567.1 1/1", 2, 25), ("@ABC1234567.4 4/1", 2, 25)],
    "SWD*EV": [("@ABC1234567.1 1/1", 2, 26), ("@ABC1234567.4 4/1", 2, 26)],
    "WD*EVR": [("@ABC1234567.1 1/1", 2, 27), ("@ABC1234567.4 4/1", 2, 27)],
    "TATLTV": [("@ABC1234567.2 2/1", 0, 0)],
    "ATLTVL": [("@ABC1234567.2 2/1", 0, 1)],
    "TLTVLA": [("@ABC1234567.2 2/1", 0, 2)],
    "LTVLAL": [("@ABC1234567.2 2/1", 0, 3)],
    "TVLALC": [("@ABC1234567.2 2/1", 0, 4)],
    "VLALCF": [("@ABC1234567.2 2/1", 0, 5)],
    "LALCFH": [("@ABC1234567.2 2/1", 0, 6)],
    "ALCFHS": [("@ABC1234567.2 2/1", 0, 7)],
    "LCFHSV": [("@ABC1234567.2 2/1", 0, 8)],
    "CFHSVV": [("@ABC1234567.2 2/1", 0, 9)],
    "FHSVVY": [("@ABC1234567.2 2/1", 0, 10)],
    "HSVVYE": [("@ABC1234567.2 2/1", 0, 11)],
    "SVVYEA": [("@ABC1234567.2 2/1", 0, 12)],
    "VVYEAI": [("@ABC1234567.2 2/1", 0, 13)],
    "VYEAIT": [("@ABC1234567.2 2/1", 0, 14)],
    "YEAITT": [("@ABC1234567.2 2/1", 0, 15)],
    "EAITTR": [("@ABC1234567.2 2/1", 0, 16)],
    "AITTRY": [("@ABC1234567.2 2/1", 0, 17)],
    "ITTRYA": [("@ABC1234567.2 2/1", 0, 18)],
    "TTRYAL": [("@ABC1234567.2 2/1", 0, 19)],
    "TRYALF": [("@ABC1234567.2 2/1", 0, 20)],
    "RYALFL": [("@ABC1234567.2 2/1", 0, 21)],
    "YALFLL": [("@ABC1234567.2 2/1", 0, 22)],
    "ALFLLR": [("@ABC1234567.2 2/1", 0, 23)],
    "LFLLRL": [("@ABC1234567.2 2/1", 0, 24)],
    "FLLRLI": [("@ABC1234567.2 2/1", 0, 25)],
    "LLRLIP": [("@ABC1234567.2 2/1", 0, 26)],
    "LRLIPK": [("@ABC1234567.2 2/1", 0, 27)],
    "PPRLPF": [("@ABC1234567.2 2/1", 1, 0)],
    "PRLPFW": [("@ABC1234567.2 2/1", 1, 1)],
    "RLPFWR": [("@ABC1234567.2 2/1", 1, 2)],
    "LPFWRY": [("@ABC1234567.2 2/1", 1, 3)],
    "PFWRYA": [("@ABC1234567.2 2/1", 1, 4)],
    "FWRYAS": [("@ABC1234567.2 2/1", 1, 5)],
    "WRYASI": [("@ABC1234567.2 2/1", 1, 6)],
    "RYASIL": [("@ABC1234567.2 2/1", 1, 7)],
    "YASILL": [("@ABC1234567.2 2/1", 1, 8)],
    "ASILLS": [("@ABC1234567.2 2/1", 1, 9)],
    "SILLST": [("@ABC1234567.2 2/1", 1, 10)],
    "ILLSTR": [("@ABC1234567.2 2/1", 1, 11)],
    "LLSTRR": [("@ABC1234567.2 2/1", 1, 12)],
    "LSTRR*": [("@ABC1234567.2 2/1", 1, 13)],
    "STRR*Q": [("@ABC1234567.2 2/1", 1, 14)],
    "TRR*QH": [("@ABC1234567.2 2/1", 1, 15)],
    "RR*QHD": [("@ABC1234567.2 2/1", 1, 16)],
    "R*QHDT": [("@ABC1234567.2 2/1", 1, 17)],
    "*QHDTL": [("@ABC1234567.2 2/1", 1, 18)],
    "QHDTLC": [("@ABC1234567.2 2/1", 1, 19)],
    "HDTLCS": [("@ABC1234567.2 2/1", 1, 20)],
    "DTLCSY": [("@ABC1234567.2 2/1", 1, 21)],
    "TLCSYS": [("@ABC1234567.2 2/1", 1, 22)],
    "LCSYSD": [("@ABC1234567.2 2/1", 1, 23)],
    "CSYSDL": [("@ABC1234567.2 2/1", 1, 24)],
    "SYSDLF": [("@ABC1234567.2 2/1", 1, 25)],
    "YSDLFR": [("@ABC1234567.2 2/1", 1, 26)],
    "SDLFRS": [("@ABC1234567.2 2/1", 1, 27)],
    "RHAYRF": [("@ABC1234567.2 2/1", 2, 0)],
    "HAYRFG": [("@ABC1234567.2 2/1", 2, 1)],
    "AYRFGA": [("@ABC1234567.2 2/1", 2, 2)],
    "YRFGAM": [("@ABC1234567.2 2/1", 2, 3)],
    "RFGAML": [("@ABC1234567.2 2/1", 2, 4)],
    "FGAMLP": [("@ABC1234567.2 2/1", 2, 5)],
    "GAMLPF": [("@ABC1234567.2 2/1", 2, 6)],
    "AMLPFC": [("@ABC1234567.2 2/1", 2, 7)],
    "MLPFCC": [("@ABC1234567.2 2/1", 2, 8)],
    "LPFCCL": [("@ABC1234567.2 2/1", 2, 9)],
    "PFCCLR": [("@ABC1234567.2 2/1", 2, 10)],
    "FCCLRG": [("@ABC1234567.2 2/1", 2, 11)],
    "CCLRGD": [("@ABC1234567.2 2/1", 2, 12)],
    "CLRGDN": [("@ABC1234567.2 2/1", 2, 13)],
    "LRGDNN": [("@ABC1234567.2 2/1", 2, 14)],
    "RGDNNT": [("@ABC1234567.2 2/1", 2, 15)],
    "GDNNTI": [("@ABC1234567.2 2/1", 2, 16)],
    "DNNTIR": [("@ABC1234567.2 2/1", 2, 17)],
    "NNTIRS": [("@ABC1234567.2 2/1", 2, 18)],
    "NTIRSV": [("@ABC1234567.2 2/1", 2, 19)],
    "TIRSVL": [("@ABC1234567.2 2/1", 2, 20)],
    "IRSVLT": [("@ABC1234567.2 2/1", 2, 21)],
    "RSVLTQ": [("@ABC1234567.2 2/1", 2, 22)],
    "SVLTQT": [("@ABC1234567.2 2/1", 2, 23)],
    "VLTQTY": [("@ABC1234567.2 2/1", 2, 24)],
    "LTQTYS": [("@ABC1234567.2 2/1", 2, 25)],
    "TQTYSE": [("@ABC1234567.2 2/1", 2, 26)],
    "QTYSEA": [("@ABC1234567.2 2/1", 2, 27)],
    "GPRDTR": [("@ABC1234567.3 3/1", 0, 0)],
    "PRDTRG": [("@ABC1234567.3 3/1", 0, 1)],
    "RDTRGS": [("@ABC1234567.3 3/1", 0, 2)],
    "DTRGST": [("@ABC1234567.3 3/1", 0, 3)],
    "TRGSTL": [("@ABC1234567.3 3/1", 0, 4)],
    "RGSTLP": [("@ABC1234567.3 3/1", 0, 5)],
    "GSTLPT": [("@ABC1234567.3 3/1", 0, 6)],
    "STLPTV": [("@ABC1234567.3 3/1", 0, 7)],
    "TLPTVL": [("@ABC1234567.3 3/1", 0, 8)],
    "LPTVLY": [("@ABC1234567.3 3/1", 0, 9)],
    "PTVLYP": [("@ABC1234567.3 3/1", 0, 10)],
    "TVLYPC": [("@ABC1234567.3 3/1", 0, 11)],
    "VLYPCF": [("@ABC1234567.3 3/1", 0, 12)],
    "LYPCFH": [("@ABC1234567.3 3/1", 0, 13)],
    "YPCFHL": [("@ABC1234567.3 3/1", 0, 14)],
    "PCFHLI": [("@ABC1234567.3 3/1", 0, 15)],
    "CFHLIA": [("@ABC1234567.3 3/1", 0, 16)],
    "FHLIAT": [("@ABC1234567.3 3/1", 0, 17)],
    "HLIATP": [("@ABC1234567.3 3/1", 0, 18)],
    "LIATPR": [("@ABC1234567.3 3/1", 0, 19)],
    "IATPRS": [("@ABC1234567.3 3/1", 0, 20)],
    "ATPRSK": [("@ABC1234567.3 3/1", 0, 21)],
    "TPRSKR": [("@ABC1234567.3 3/1", 0, 22)],
    "PRSKRT": [("@ABC1234567.3 3/1", 0, 23)],
    "RSKRTH": [("@ABC1234567.3 3/1", 0, 24)],
    "SKRTHP": [("@ABC1234567.3 3/1", 0, 25)],
    "KRTHPL": [("@ABC1234567.3 3/1", 0, 26)],
    "RTHPLL": [("@ABC1234567.3 3/1", 0, 27)],
    "ALEIRA": [("@ABC1234567.3 3/1", 1, 0)],
    "LEIRAG": [("@ABC1234567.3 3/1", 1, 1)],
    "EIRAGV": [("@ABC1234567.3 3/1", 1, 2)],
    "IRAGVR": [("@ABC1234567.3 3/1", 1, 3)],
    "RAGVRS": [("@ABC1234567.3 3/1", 1, 4)],
    "AGVRSQ": [("@ABC1234567.3 3/1", 1, 5)],
    "GVRSQL": [("@ABC1234567.3 3/1", 1, 6)],
    "VRSQLF": [("@ABC1234567.3 3/1", 1, 7)],
    "RSQLFY": [("@ABC1234567.3 3/1", 1, 8)],
    "SQLFYT": [("@ABC1234567.3 3/1", 1, 9)],
    "QLFYTP": [("@ABC1234567.3 3/1", 1, 10)],
    "LFYTPV": [("@ABC1234567.3 3/1", 1, 11)],
    "FYTPVF": [("@ABC1234567.3 3/1", 1, 12)],
    "YTPVFI": [("@ABC1234567.3 3/1", 1, 13)],
    "TPVFIL": [("@ABC1234567.3 3/1", 1, 14)],
    "PVFIL*": [("@ABC1234567.3 3/1", 1, 15)],
    "VFIL*Q": [("@ABC1234567.3 3/1", 1, 16)],
    "FIL*QH": [("@ABC1234567.3 3/1", 1, 17)],
    "IL*QHP": [("@ABC1234567.3 3/1", 1, 18)],
    "L*QHPG": [("@ABC1234567.3 3/1", 1, 19)],
    "*QHPGV": [("@ABC1234567.3 3/1", 1, 20)],
    "QHPGVS": [("@ABC1234567.3 3/1", 1, 21)],
    "HPGVSA": [("@ABC1234567.3 3/1", 1, 22)],
    "PGVSAL": [("@ABC1234567.3 3/1", 1, 23)],
    "GVSALI": [("@ABC1234567.3 3/1", 1, 24)],
    "VSALIL": [("@ABC1234567.3 3/1", 1, 25)],
    "SALILS": [("@ABC1234567.3 3/1", 1, 26)],
    "ALILSY": [("@ABC1234567.3 3/1", 1, 27)],
    "PSRYAR": [("@ABC1234567.3 3/1", 2, 0)],
    "SRYARE": [("@ABC1234567.3 3/1", 2, 1)],
    "RYAREY": [("@ABC1234567.3 3/1", 2, 2)],
    "YAREYA": [("@ABC1234567.3 3/1", 2, 3)],
    "AREYAP": [("@ABC1234567.3 3/1", 2, 4)],
    "REYAPN": [("@ABC1234567.3 3/1", 2, 5)],
    "EYAPNC": [("@ABC1234567.3 3/1", 2, 6)],
    "YAPNCF": [("@ABC1234567.3 3/1", 2, 7)],
    "APNCFI": [("@ABC1234567.3 3/1", 2, 8)],
    "PNCFIP": [("@ABC1234567.3 3/1", 2, 9)],
    "NCFIPL": [("@ABC1234567.3 3/1", 2, 10)],
    "CFIPLF": [("@ABC1234567.3 3/1", 2, 11)],
    "FIPLFS": [("@ABC1234567.3 3/1", 2, 12)],
    "IPLFSS": [("@ABC1234567.3 3/1", 2, 13)],
    "PLFSSY": [("@ABC1234567.3 3/1", 2, 14)],
    "LFSSYS": [("@ABC1234567.3 3/1", 2, 15)],
    "FSSYSN": [("@ABC1234567.3 3/1", 2, 16)],
    "SSYSNT": [("@ABC1234567.3 3/1", 2, 17)],
    "SYSNTP": [("@ABC1234567.3 3/1", 2, 18)],
    "YSNTPE": [("@ABC1234567.3 3/1", 2, 19)],
    "SNTPE*": [("@ABC1234567.3 3/1", 2, 20)],
    "NTPE*A": [("@ABC1234567.3 3/1", 2, 21)],
    "TPE*AH": [("@ABC1234567.3 3/1", 2, 22)],
    "PE*AHS": [("@ABC1234567.3 3/1", 2, 23)],
    "E*AHSS": [("@ABC1234567.3 3/1", 2, 24)],
    "*AHSSS": [("@ABC1234567.3 3/1", 2, 25)],
    "AHSSSP": [("@ABC1234567.3 3/1", 2, 26)],
    "HSSSPI": [("@ABC1234567.3 3/1", 2, 27)],
    "TATHPT": [("@ABC1234567.5 5/1", 0, 0)],
    "ATHPTL": [("@ABC1234567.5 5/1", 0, 1)],
    "THPTL*": [("@ABC1234567.5 5/1", 0, 2)],
    "HPTL*E": [("@ABC1234567.5 5/1", 0, 3)],
    "PTL*ED": [("@ABC1234567.5 5/1", 0, 4)],
    "TL*EDI": [("@ABC1234567.5 5/1", 0, 5)],
    "L*EDIN": [("@ABC1234567.5 5/1", 0, 6)],
    "*EDING": [("@ABC1234567.5 5/1", 0, 7)],
    "EDINGD": [("@ABC1234567.5 5/1", 0, 8)],
    "DINGDR": [("@ABC1234567.5 5/1", 0, 9)],
    "INGDRC": [("@ABC1234567.5 5/1", 0, 10)],
    "NGDRCT": [("@ABC1234567.5 5/1", 0, 11)],
    "GDRCTN": [("@ABC1234567.5 5/1", 0, 12)],
    "DRCTNR": [("@ABC1234567.5 5/1", 0, 13)],
    "RCTNRA": [("@ABC1234567.5 5/1", 0, 14)],
    "CTNRAD": [("@ABC1234567.5 5/1", 0, 15)],
    "TNRADA": [("@ABC1234567.5 5/1", 0, 16)],
    "NRADAH": [("@ABC1234567.5 5/1", 0, 17)],
    "RADAHY": [("@ABC1234567.5 5/1", 0, 18)],
    "ADAHYF": [("@ABC1234567.5 5/1", 0, 19)],
    "DAHYFT": [("@ABC1234567.5 5/1", 0, 20)],
    "AHYFT*": [("@ABC1234567.5 5/1", 0, 21)],
    "HYFT*V": [("@ABC1234567.5 5/1", 0, 22)],
    "YFT*VV": [("@ABC1234567.5 5/1", 0, 23)],
    "FT*VVG": [("@ABC1234567.5 5/1", 0, 24)],
    "T*VVGG": [("@ABC1234567.5 5/1", 0, 25)],
    "*VVGGS": [("@ABC1234567.5 5/1", 0, 26)],
    "VVGGSR": [("@ABC1234567.5 5/1", 0, 27)],
    "PPRILP": [("@ABC1234567.5 5/1", 1, 0)],
    "PRILPC": [("@ABC1234567.5 5/1", 1, 1)],
    "RILPCK": [("@ABC1234567.5 5/1", 1, 2)],
    "ILPCKR": [("@ABC1234567.5 5/1", 1, 3)],
    "LPCKRI": [("@ABC1234567.5 5/1", 1, 4)],
    "PCKRIS": [("@ABC1234567.5 5/1", 1, 5)],
    "CKRISM": [("@ABC1234567.5 5/1", 1, 6)],
    "KRISMA": [("@ABC1234567.5 5/1", 1, 7)],
    "RISMAI": [("@ABC1234567.5 5/1", 1, 8)],
    "ISMAIG": [("@ABC1234567.5 5/1", 1, 9)],
    "SMAIGV": [("@ABC1234567.5 5/1", 1, 10)],
    "MAIGVQ": [("@ABC1234567.5 5/1", 1, 11)],
    "AIGVQT": [("@ABC1234567.5 5/1", 1, 12)],
    "IGVQTE": [("@ABC1234567.5 5/1", 1, 13)],
    "GVQTEL": [("@ABC1234567.5 5/1", 1, 14)],
    "VQTELM": [("@ABC1234567.5 5/1", 1, 15)],
    "QTELMP": [("@ABC1234567.5 5/1", 1, 16)],
    "TELMPT": [("@ABC1234567.5 5/1", 1, 17)],
    "ELMPTI": [("@ABC1234567.5 5/1", 1, 18)],
    "LMPTIS": [("@ABC1234567.5 5/1", 1, 19)],
    "MPTISR": [("@ABC1234567.5 5/1", 1, 20)],
    "PTISRK": [("@ABC1234567.5 5/1", 1, 21)],
    "TISRK*": [("@ABC1234567.5 5/1", 1, 22)],
    "ISRK*W": [("@ABC1234567.5 5/1", 1, 23)],
    "SRK*WE": [("@ABC1234567.5 5/1", 1, 24)],
    "RK*WEG": [("@ABC1234567.5 5/1", 1, 25)],
    "K*WEGR": [("@ABC1234567.5 5/1", 1, 26)],
    "*WEGRV": [("@ABC1234567.5 5/1", 1, 27)],
    "RHASYL": [("@ABC1234567.5 5/1", 2, 0)],
    "HASYLV": [("@ABC1234567.5 5/1", 2, 1)],
    "ASYLVR": [("@ABC1234567.5 5/1", 2, 2)],
    "SYLVRG": [("@ABC1234567.5 5/1", 2, 3)],
    "YLVRGY": [("@ABC1234567.5 5/1", 2, 4)],
    "LVRGYQ": [("@ABC1234567.5 5/1", 2, 5)],
    "VRGYQW": [("@ABC1234567.5 5/1", 2, 6)],
    "RGYQWR": [("@ABC1234567.5 5/1", 2, 7)],
    "GYQWRS": [("@ABC1234567.5 5/1", 2, 8)],
    "YQWRSV": [("@ABC1234567.5 5/1", 2, 9)],
    "QWRSVY": [("@ABC1234567.5 5/1", 2, 10)],
    "WRSVYK": [("@ABC1234567.5 5/1", 2, 11)],
    "RSVYKQ": [("@ABC1234567.5 5/1", 2, 12)],
    "SVYKQS": [("@ABC1234567.5 5/1", 2, 13)],
    "VYKQS*": [("@ABC1234567.5 5/1", 2, 14)],
    "YKQS*C": [("@ABC1234567.5 5/1", 2, 15)],
    "KQS*CP": [("@ABC1234567.5 5/1", 2, 16)],
    "QS*CPL": [("@ABC1234567.5 5/1", 2, 17)],
    "S*CPLF": [("@ABC1234567.5 5/1", 2, 18)],
    "*CPLFH": [("@ABC1234567.5 5/1", 2, 19)],
    "CPLFHV": [("@ABC1234567.5 5/1", 2, 20)],
    "PLFHVS": [("@ABC1234567.5 5/1", 2, 21)],
    "LFHVSS": [("@ABC1234567.5 5/1", 2, 22)],
    "FHVSSG": [("@ABC1234567.5 5/1", 2, 23)],
    "HVSSGR": [("@ABC1234567.5 5/1", 2, 24)],
    "VSSGRV": [("@ABC1234567.5 5/1", 2, 25)],
    "SSGRVA": [("@ABC1234567.5 5/1", 2, 26)],
    "SGRVAC": [("@ABC1234567.5 5/1", 2, 27)],
}
