(* Content-type: application/vnd.wolfram.cdf.text *)

(*** Wolfram CDF File ***)
(* http://www.wolfram.com/cdf *)

(* CreatedBy='Mathematica 9.0' *)

(*************************************************************************)
(*                                                                       *)
(*  The Mathematica License under which this file was created prohibits  *)
(*  restricting third parties in receipt of this file from republishing  *)
(*  or redistributing it by any means, including but not limited to      *)
(*  rights management or terms of use, without the express consent of    *)
(*  Wolfram Research, Inc. For additional information concerning CDF     *)
(*  licensing and redistribution see:                                    *)
(*                                                                       *)
(*        www.wolfram.com/cdf/adopting-cdf/licensing-options.html        *)
(*                                                                       *)
(*************************************************************************)

(*CacheID: 234*)
(* Internal cache information:
NotebookFileLineBreakTest
NotebookFileLineBreakTest
NotebookDataPosition[      1063,         20]
NotebookDataLength[    468892,       9629]
NotebookOptionsPosition[    459433,       9369]
NotebookOutlinePosition[    459954,       9392]
CellTagsIndexPosition[    459911,       9389]
WindowFrame->Normal*)

(* Beginning of Notebook Content *)
Notebook[{

Cell[CellGroupData[{
Cell["\<\
Elemento isoparam\[EAcute]trico de 4 Nodos:\
\>", "Section", "PluginEmbeddedContent"],

Cell[CellGroupData[{

Cell["\<\
Ing. Jorge De La Cruz. Msc.
An\[AAcute]lsis por medio de Elementos Finitos.\
\>", "Subsection", "PluginEmbeddedContent"],

Cell["Definiendo las coordenadas nodales:", "Text", "PluginEmbeddedContent"],

Cell[BoxData[{
 RowBox[{
  RowBox[{"x1", "=", "0.0"}], ";", " ", 
  RowBox[{"y1", "=", "0.0"}], ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{"x2", "=", "5.0"}], ";", " ", 
  RowBox[{"y2", "=", "0.0"}], ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{"x3", "=", "5.0"}], ";", " ", 
  RowBox[{"y3", "=", "5.0"}], ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{"x4", "=", "0.0"}], ";", 
  RowBox[{"y4", "=", "5.0"}], ";"}]}], "Input", "PluginEmbeddedContent"],

Cell["Graficando el elemento:", "Text", "PluginEmbeddedContent"],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{"ListLinePlot", "[", 
  RowBox[{
   RowBox[{"{", 
    RowBox[{
     RowBox[{"{", 
      RowBox[{"x1", ",", "y1"}], "}"}], ",", 
     RowBox[{"{", 
      RowBox[{"x2", ",", "y2"}], "}"}], ",", 
     RowBox[{"{", 
      RowBox[{"x3", ",", "y3"}], "}"}], ",", 
     RowBox[{"{", 
      RowBox[{"x4", ",", "y4"}], "}"}], ",", 
     RowBox[{"{", 
      RowBox[{"x1", ",", "y1"}], "}"}]}], "}"}], ",", 
   RowBox[{"Filling", "\[Rule]", "Bottom"}]}], "]"}]], "Input", \
"PluginEmbeddedContent"],

Cell[BoxData[
 GraphicsBox[{{}, 
   GraphicsComplexBox[{{0., 0.}, {5., 0.}, {5., 5.}, {0., 5.}, {0., 
    0.}}, {{{}, 
      {Hue[0.67, 0.6, 0.6], Opacity[0.2], EdgeForm[None], 
       GraphicsGroupBox[PolygonBox[{{2, 5, 4, 3}}]]}, {}, {}}, {{}, {}, 
      {RGBColor[0.24720000000000014`, 0.24, 0.6], 
       LineBox[{1, 2, 3, 4, 5}]}}}], {}},
  AspectRatio->NCache[GoldenRatio^(-1), 0.6180339887498948],
  Axes->True,
  AxesLabel->{None, None},
  AxesOrigin->{0, 0},
  Method->{},
  PlotRange->{{0, 5.}, {0, 5.}},
  PlotRangeClipping->True,
  PlotRangePadding->{{0.1, 0.1}, {0.1, 0.1}}]], "Output", \
"PluginEmbeddedContent"],

Cell[BoxData[
 GraphicsBox[{{}, 
   GraphicsComplexBox[{{0., 0.}, {4., 0.}, {4., 1.}, {1., 1.}, {0., 
    0.}}, {{{}, 
      {Hue[0.67, 0.6, 0.6], Opacity[0.2], EdgeForm[None], 
       GraphicsGroupBox[PolygonBox[{{2, 5, 4, 3}}]]}, {}, {}}, {{}, {}, 
      {RGBColor[0.24720000000000014`, 0.24, 0.6], 
       LineBox[{1, 2, 3, 4, 5}]}}}], {}},
  AspectRatio->NCache[GoldenRatio^(-1), 0.6180339887498948],
  Axes->True,
  AxesLabel->{None, None},
  AxesOrigin->{0, 0},
  Method->{},
  PlotRange->{{0, 4.}, {0, 1.}},
  PlotRangeClipping->True,
  PlotRangePadding->{{0.08, 0.08}, {0.02, 0.02}}]], "Output", \
"PluginEmbeddedContent",
 GeneratedCell->False,
 CellAutoOverwrite->False]
}, Open  ]],

Cell["\<\
Definiendo las funciones de forma en el espacio (\[Xi],\[Eta]):\
\>", "Text", "PluginEmbeddedContent"],

Cell[BoxData[{
 RowBox[{
  RowBox[{"N1", "=", 
   RowBox[{
    RowBox[{"1", "/", "4"}], 
    RowBox[{"(", 
     RowBox[{"1", "-", "\[Xi]"}], ")"}], 
    RowBox[{"(", 
     RowBox[{"1", "-", "\[Eta]"}], ")"}]}]}], ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{"N2", "=", 
   RowBox[{
    RowBox[{"1", "/", "4"}], 
    RowBox[{"(", 
     RowBox[{"1", "+", "\[Xi]"}], ")"}], 
    RowBox[{"(", 
     RowBox[{"1", "-", "\[Eta]"}], ")"}]}]}], ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{"N3", "=", 
   RowBox[{
    RowBox[{"1", "/", "4"}], 
    RowBox[{"(", 
     RowBox[{"1", "+", "\[Xi]"}], ")"}], 
    RowBox[{"(", 
     RowBox[{"1", "+", "\[Eta]"}], ")"}]}]}], ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{"N4", "=", 
   RowBox[{
    RowBox[{"1", "/", "4"}], 
    RowBox[{"(", 
     RowBox[{"1", "-", "\[Xi]"}], ")"}], 
    RowBox[{"(", 
     RowBox[{"1", "+", "\[Eta]"}], ")"}]}]}], ";"}]}], "Input", \
"PluginEmbeddedContent"],

Cell["Graficando las funciones de forma:", "Text", "PluginEmbeddedContent"],

Cell[CellGroupData[{

Cell[BoxData[{
 StyleBox[
  RowBox[{"Plot3D", "[", 
   RowBox[{"N1", ",", 
    RowBox[{"{", 
     RowBox[{"\[Xi]", ",", 
      RowBox[{"-", "1"}], ",", "1"}], "}"}], ",", 
    RowBox[{"{", 
     RowBox[{"\[Eta]", ",", 
      RowBox[{"-", "1"}], ",", "1"}], "}"}]}], "]"}], 
  "Input"], "\[IndentingNewLine]", 
 RowBox[{"Plot3D", "[", 
  RowBox[{
   RowBox[{"{", 
    RowBox[{"N1", ",", "N2", ",", "N3", ",", "N4"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"\[Xi]", ",", 
     RowBox[{"-", "1"}], ",", "1"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"\[Eta]", ",", 
     RowBox[{"-", "1"}], ",", "1"}], "}"}]}], "]"}]}], "Input", \
"PluginEmbeddedContent"],

Cell[BoxData[
 Graphics3DBox[GraphicsComplex3DBox[CompressedData["
1:eJx1nAucjeX2x0dMRe4hhpAGcRw5JTU6zUMo2y0lx52KisNGU2kSqXEZ1Bza
ZMb1ROWSS26N2/CajcmdmWFcxmVcBsOQcifNf57Zz3e99vN3+pzPZ3/mdx7v
Xu9vr3c961nrt94n3hn0xrsPhISEVAkNCSmc//nE0bYr8vJ+c/hsF/hU7320
dGqi94LgoW//s05m4kX10yOPfuHNPCP42tN1/5u/Tl2a8/F74Z6Tgg/uW6m8
LzxHvVtwnUzBa+U+9FX+vzfr0wQ/Ouh6iMd32uB+wX1XTg/Jv64KMf+BT7w0
91Be3nHWK/A710Zuyf8+ri9437tvLfNmZmCP4PuKRM7U64z9gn9RsVHZVd4d
3K/g/6p9OUd/n+FH8G4D3hxZ07OcvwXv5Gz5bmy9DyP5G17hk7/5nJ6ztFdM
/AXhGbxBVM0pdaJzhGfwLben7tzb+azwDN51ZMnC0RHZwjP45UdGRlQLOyU8
g4+efGNQyu0s4Rn82cYH22QmHhWewXu0Tl5TJ/qQ8Aw+tvtPtaIj9gvP4Cu8
kyal3N4rPIOvvXszpG70NuEZvGT1HYdTbm8UnsE3jWrc79fbS4Vn8DVlLk3J
3fWB8Ax/+CfrwOGTv/n8Ia9a3IyYs8IzeMW4sFGl4rKFZ/C4sPKfxcSfEp7B
H5hfKura7BPCM3juv8ul7e18XHgG3z+66EbXb/1i14bv/lo8I+aA8Aw+b92V
6THx6cIzeI3FqY+mdt4jPINfObnutxkxvwrP4OP+0SY+tbMjPIP/q3HTWWmd
fxaewcP6VIhrP2iw8AxPPO/wDI5/8u/B4ZO/+RxQbfqL7ZJPCc/gG5dVvvL1
tRPCM3i55jN+2lEnS3gG37soYnO5pCPCM/gnqY2f9YUfEp7Bq15/cU6puP3C
M/ivYS+VibuWKjyDH/18WEj5pJ3CM3i9rZ8fKxW3RXgGv17jgU/LJyUJz+C1
+5f4d4WkRcKzrB8w/MF+TwwUnuGD+AnP4Dzv8Awe7J8nBYdP/ubzSssa/eM7
ZAnP4G03vp7SrLfEAcFXXOqZ2LnaYeEZvNLjA+bqeAvP4CNafzpF+y08g5/5
dMyY+A57hWfwqe1rX+xcbbvwLPYubrorJn6T8Azeu9zqrl2qrRWewW9MP/Rc
12o/Cc/gq8vs/qNBn/6WP6fJfhTsz2kSP4P9Oc163k8Kjn8G+3OaxWea4B/e
/Pfzl6MOW/6c5sQ3GDVfx4Fgf05zkvrOrOTx7bP8Oc058V3i+GuzUy1/TnOK
L//vjWa9d1n+nOb8q0baIR0fgv05zfmvr0jStdnJlj+nOamhA/85M2aV5c9p
ztadPxaaFTPP8uc057PEqsdvzu1r8eyX/T2YZ7/sR8E8+yV+BvPsl+c9mGe/
+Gcwz37hM5hnv/NZ1IJpWdn7LJ79Ts+Vu792efTLvxuZe/buXom38v87rQY/
c1L7bTDPfqfs1WEpWdmbLZ79zrLf/tNRx+Fgnv2OuvFk1RPZKy2e/U7C3Q92
ncj+weLZ70R2G7RjY867wjP3Qb4Ez+Ds7/AMHrwfnRSc+AnP4MHPe5rgln8K
Dp8h9/1vr/AM0jz+nTOJ+fkYPIOPq79lpzczRXgG372l9op74rvg4WvLPbvK
u054Bl907lBeeH7+Bs/gfetHTq3pmSM8uzZuXDW2Xm/bnxX5p+XPinzJ8mcV
vL+LPyv2I8uflRU/BQ9+3sWfFf5p+bMK5lP8WYW+fHS3Gwdcf9+a22mVjreW
P6tyT/+otN9a/qziR00rreOD5c/q2ZXVNqXm5xWWP6vBG8ptT+v8X9uflf+H
0j+0H/SWHZ8V+bwVnxX5pxWfFfmSFZ+Vtb8Lzn5kxWdF/LTis+J5t+Kzwj+t
+Cx8WvFZPbTw/Oz4Dhvt+Kzi0m5X1PuaFZ9VsdeqX4qJX2HHZ1Uy/qsxXfL3
Oys+qzLVF7/dtdo0Oz6rYpP81f/Rp5udbyjOR1a+ocjnrXxDWfmn4ORLVr6h
gvd3yTdU8H4k+YYifsIzePDzLvmG5Z+SbwifVr6hQoe0/yMrP95a+Ya6/YzO
k3+28w11/tMbHh2HrXxDrfiyxa0T2ZPtfENNPNL71sacjnb+rDhvWvmzss5H
gpPPW/mzIv+08mdFvmTlz4r93cqfFftRsD+fVMRPK3+W593Kn8U/rfxZ+LTy
Z3X8105vzMzPH6z8WY1bkxKm44OVP6u5KeueSO88wc6fVflaoWntB7W3z4PK
Or8LznnTOg8qzkfwDE4+b50HFfknPIOTL8EzOPs7PIOzH8EzOPETnsV+87xb
50HxT3iR+zJ8WudBFTli1mztt9Z5UJ0Zl9v+ZPYo+zyonj3epmNyTgu7vqGo
h1j1DcX53apvKM6b8AzO+ciqbyjyeXgGJ/+EZ3DyJXgGZ3+HZ3D2I3gGJ35a
9Q153uEZHP+06hvCp1XfUG2nNQ45lf2RXd9Q0aerlvDnPCc88/9TX7LqdYp6
iFWvU5zf4Rmc86ZVr1Ocj+AZnHwensHJP+HZzZsC+RI8g7O/wzM4+5FVr5P4
Cc/gPO9WvU7806rX2XwK/lfm+flNbl//5+anRtTKTLzo8PnSyR6V2iVflLoe
eN1qLQZ6My8p8zefjqdG5Vk3Zs+W64KPPXdoeuPhUZHW9ZW1XvDe16sk1sr3
T+v66uQfUeuzs9+Uv/meQr/5h23KqSnX53q2neDW9eV3xE7Lfrn+HzFnBuRf
z+EzekTbusWOXVLmbz6dp5qsHvXNtW+UtV5ZuHyWrH/6+QvZDcU+/t0/P982
MdH7m6zj/x9VbEnlOh6fXB97f0t/bdaJfa798r1FHumeklM+0lov17fWO9Z1
xH7re+Xfcf3hT5xOystz67pjO5+M9mZeUIfeHZqSl3dWvndXs5+T89cJzvpH
2o/dmpeXo577qfTevLxTsr74jeNBOOvjXq1f0Rd+TtUoWHdc1rcuWHdWcNaX
VvveD/ecUYcLvveArP+qR9MDeXnZgrN+UqOhqxK9p1V4wbo9sn5H6ahj+XYI
zvoK9as/5PGdVLULcEfWF908Jzsv74SJwwlJ4N+Y8ynruU5CzZR/ZSZmyXrw
hubcWjmwXvyhZMeHz+XlZQrO+hfqhOfm5R1SVQN2yvq4df2u59+n4KyPrDj4
q3BPhjoQ4EHWF6+xo0hIyD7BWd/8oaS6Hl+aqhbgWdZ/Nbbeo/o8C856X26p
yZmJu1TDwO8o68scq9wsJGSr4Kyftn1eWrgnRR0M+Il7vwuP9cqPrYKzvtXr
I2ocSVyvhgX8UNY/dvnlWyEhSwVn/bqwnupo4mKF37JvfFn9leczE3PEP8Ff
7f72Hv19+CF4iYRhBX6Fv4Gnp8fntfSdFr8Cn1pqRbw386T4D3iv1ruf9oWf
EH8Ab1UyfYI385j4A3iGqVfwu4O/c1ddvvf3Bb/40ZDqui7B7wg+NHdRe194
qvxe4MW3phTS+Tm/i9TbVz37jse3WfgH/7hN09E6T4ZncE/HRl+v9i4UntmH
h55bcK2l76zwDD6yb+MM/TzCM3jcue2J+nmEZ/ApfbsV8AnP4N+du/BJuOeE
8AzeIfXm0kSv+9xJ/pVaKSrcc0R4Bo829R94Bh+f2u1qXt5+4Rl8ZuqwX7yZ
acIz+Fvfv3gx0btHeAYf8/3wQb7wbcKz5N3fO797M/3Cs9zXq4XP+cLXCM/g
ByeUPDQpfIHwTF7z4dISKxO92cKz1MM/evSo/j54Bi8RUelBHcfgWertd6sW
+Cc8S518pKeAT3gGP/tIxwc8viPCM/gXk996MzPxkPAMHmbqafAs9fYfP7mZ
l5cuPMv31h/ZSvstPIPv8j4e4fHtFJ6lbjei4VFdl4Bn8Gcmtv7SF54sPIPv
fHlLq4GZicKz5OPbz9QclDnX8uc0Z2D0B/n+c8ry5zTnry9yT+n4H+zPaU7c
uPc76jgf7M9pzovhW+L08x7sz2lO9+fWDdb+GezPac7wV5Z10HwG+3OaM6PT
vEa+8AOWP7v182B/TnOOfjrpz7y8VMuf05zb5R91vJm7LX/Ot/OFvFGZidss
f05zhnW90ErHh2B/TnOG3FryxMDMDZY/pzk3X4pd7PGttPw5zblaxIlq5fvB
4jn/fNNkfwPth8E8+52Pwr7ar+0I5tnvVF7Z36P3zWCe/Y6/bZv1On4G8+x3
+p2t9w/9vAfz7HdKf1niR+2fwTz7ndVhlypqPoN5duvnwTz7nQfb/VxQ7w3m
2e902VMtR9cng3n2O/MfOLFX+20wz37n+nNzVutzRTDPfic2LaWJjsPBPPud
8AhPaX0+DebZ78SoeEefB+GZc8yGgjznqPAMnvvMS6d0vgHP4JUL8pBDwjN4
y0/Sc/X9wDP4J1Nn/a7jJzyD/1iQn6QLz+DpRxve0f4JzyFB/+0VnkG6N59S
NCRkp/AMvrogP9kqPIOX3+WvEhKyWXgG/6BMl1o634Bn8ISCPGS18AzevuEj
H96bb4CXHjqoeKGQ72x/Vsu+i+imn1/Ln9V7rX8r2Mctf1Zh138Yo/M0y5/V
nu+6VtF+aPmzGtW69HK9H1n+rBpf3/Kqjp+WP6sjTW+s1c+75c+qhamfW/6s
luTuKqP5tPxZVXw54iUdByx/VjHxP/TV8dbyZ3XgULWntN9a/qwWjLr+uY4P
lj+rmz2nLVyVn1dY/qz+2W/C8tXe6XZ8Vj80KHQw0bvfjs/qeO/Zf2l/s+Kz
qhzftKbOb634rDptP9Fa82bFZ1XzTuFsvb9b8VnNeLv5M3o/suKzKrt11Agd
P634rMab+rkVn1Whb0PDtH9a8VlF32nxvubTis/q2Kbi/Ty+JDs+q10T4wrr
fc2Kz2ptj5KzfOHL7Pis2nes+sHA/P3Ois+q8zHVZlDmFDvfUPEbt76p+bHy
DdV1WIfi+jm18g31/OCkWH0usPINte5s6nXtb1a+oZr0OvuuzpesfEOlZPy5
T/Nm5RuqVbuyzbUfWvmG1M+tfEN1eOmlGvp5t/IN9XO/H2to/7TyDbWm894Q
zaeVbyj/q3eO6XqmlW+ob71thmq/tfINNb/1nEY6Dlv5hhrS4XJ2viPa+bN6
dHb177X/WPmz+rHOIw3182jlz+r55dc26/OUlT+rbY2zOmp+rPxZddu0/YzO
P638WV1q/csn2t+s/FmFLFBl9P4e7M8nVU1TP7fyZ+XZkfOmjp9W/qwGHtkw
RD/vVv6sfBcnJWj/tPJnNfBu6keaTyt/VonVhq325OcPVv6s8ka2ubkqPz5Y
+bN68pumeau9X9nnQXVj488z9P1a50EV33XQa9pPrPOg8p7rUUSfQ63zoFpc
+Wg/zYN1HlQX23Xfo/N56zyo/h6T2VD7lXUeVN5fuk7T/MCzXT+3zoMqt3KX
d/V+ZJ0HVedZS7/VvFnnQXWwyM9v6efdOg+qzv0X/037p3UeVEu+73RR82md
B9WSwboe7rPPg6rrZ3OeqeWJUVa9SCUmD683Kf/3tepCqsxj8/Zr+636jzp0
8Lu7R/LP71adR00YeeU1fd606jnqlfqvzNHnI6tuo/48mHBV369Vh1HLR154
Reef8GzXz616i6p6aOKFkPz93aqrqAUP/fqZ3o+s+ol6rMeA6jp+WnUSNWp5
mS2aN6seoo6XGpiu/dOqe6jVa1dX0Xxa9Q0VvuPBnbU8UcIz9TrvQ3Uqbrz1
QaRVl1Mxr3094UbK4Eir/qbqP328RM+nBkVadTbV2fnxz53+AZFWPS3/OgMu
vNjj35FW3UwtOv7M4Z9uvB9p1cfUgUG3tlbyUT+X+pjUz616l6o7cczcGylv
RVp1LbWgybS6u/zdI636lcp7q/bDC290irTqVOrNL1dmj63XIdKqR6l2zUaO
2eVvE2nVndTNkTObjqvXTOqcfG6qXaXdtSdyxW+/3f1YVktfrpP99wu398ee
F79d/Fnrwr7wi86Xb0y8MzfLrVtShz5z7dHioU+fFb89sn+5J9xzySlbcUi9
y1HZ4rdXTJ324TWJ3fsuOiV+W6xBmK7nOpsP7p/hzXTrje//j/rh58b+zXnr
h4d7DovfYv8Z/1u7Y1PceiD2J4zcGBn6dLr4LfavGLl2ZkToHvHbTGP/v+54
7+q+LXxif6MjjSL6LXLEP5NXv/J9ovesU8zUgeGTzxqNNxfKiD0nfP5edklO
S98554Hvj7WKCHXrvUU3nijuC89xfhp/4Ps60W5dt7q3fIPMxBwn7Pclxedl
ufXbF8I8HcI95503MvqMXhbp1m9f+3XYEG9mPt74VgmPz63HLvwf9Te/sf9S
+qEvwz0Hhc/D+EW193svKrVf+MT+3X2X9WmX7NbfsL9+rw8e0foB+MT+42dv
fZ28Y4vwif19Dtyc3ry3W/9sb+yf/NOpiHvrnKHtF4zzZp50ClfMHePNdOtv
j989sijRe8rp+crMhzNizwjPfMZcqTdtR51s4bltp+ZXWvpOO16jX4XnPkWi
K/jCs53q28/m9mjn1jmHLVsYkZmY7cy4UnbJ7b1ZwvOknse7h3vOOBuej/ir
pc+tc/72P+pvDxr7K8y/Nlqfg+C5qrH/lWHem3s7pwvP2F+k3ZzxWt8Cz9i/
oPFTn16O2i48Y//7+0ZXOVPOrXNi/4S0tX9VC3PrnJON/TWy1425t845qeA5
OuL8efy1Fvn/TnjeVPC8HHWmvzGkVEbsabEv8Fwccx5o/llOj3anhGc+G6zs
mR5VzK1zzp9++M+Wviyn9vfVd5VLyhKeN7YsWd0XfsLZ+UTM4tiUY8LzgWtN
m2UmnnAO3x6Rpc/vwfva/6+/fWvsV8PKjA/37BOeNxv7rx38ME7rr+AZ+4e8
V+vt5B27hGe5jwrHz+yo49b5sf+3jWUWFevp1jmxf92uYW/f2btaeD5o7F9a
bXL6vXXOKgW/e7pzfd6Kwh6fW39rU/D77nMO9nytQkbsSeE58Dvud2Imb3vh
ctQJsS/wXGc4V40eGJ757H5zedf4DseE58Bzesh5u/uxvypdyBSeA8/jYefa
+hR/otetc142da31Vv2tqrG/X/u//0efq+C5rbG/Vf+rzRtm7BGesf90g4cX
d7i8XXjG/lm/f/vmvKwtwvMhY/8ni2e++PQ8t2+C/We2Rm7OiHXrnBHG/siu
Oj9x65ytCuzZ7Fwe9Fj+fZ2w8uEUp1z/px7PiM0SnosVXH+bU2T706VDn3b9
oG+Bv+1wBq3bERodcVR4Tinwq13Ov38Krxz6dKbwzOf1KY3bLyrl1jmTC/wk
1Slj6nLw/Lupa/Wy6m/tjP2+cE/+/9z62yJj/4LwQoOG990pPGP/63dybmaX
2yo8Y3+lXW38Wl8Ez1uM/SENO0yZl7VeeMb+RzeH1F4e6dY5/cb+XZObe++t
c2YV2O9z9hbosY9Z57s5zpxLD7RcVOqo8DywwJ6FTpH+O5fUic4Unq8WrF/u
XD/9xnPJOw4Jz0sL7jfR+aT7lMyoYgfFvgEF11nnjE+dMjuyoVvn5LNe8zdj
dRwIrr8lJBWy6m9Fvg/Y/9Hv2+rpugc83zT2P3ysxNSU2279bbCxf8OWHc2j
I7YIz9g/aX7hBj3bJQvPS4z9pUv0/+Vy1Frh2WvsV45/XNy1ZcJzPWN/tz1V
1t5b5ww39m4y+U9wvSL/vPrQLx21v8Hzj4HvVYWL92vbLvmA8HwuYKcaVvJw
n2I9M4TnzoH7UqNLPjxD66vhOfD7zlFVi1/5Q+9T8Lsv8LurVx+a3i8z0a1z
XjB1rVes+hufKVMLPefxufW3t4z9fePDNmldHDz/YOz/28TFE3W8hWfsr/Lw
xy9+WGyD8NzD2P/EmZh9EaGrhOdjxv59HxU9pHVZ4s/G/pyi1x68t84ZiOep
KnHC6Nu63gvP1Y29k01eBM8bA8+RKrI4NDeq2D7h+Z3Ac6cKv/D25A6X04Tn
IoHnVPk2ftpd51HwPDfwXKsVr7zWvPfwvcJzy0AcUFsmb4oP97h1ztz/UX/b
YOyvZOqZkm+Yz60Fek6/8Iz9DX6puiYi1BGesb+Mc7JC3eg1wnOosf8ZZ9oV
rRuEZ+yv1z9u1tiUn4TnV4z9F8dOm3dvnfP5QNxWC6q1/UT7Dzw/FIjzqtQL
lT5YVCpVeMZe8iJ4nh/YR1T5TXtjdT4Pz58G9h3VNDsjIjpil/Ac2Bf2qQkP
dPrbhWY7hOfKgX1N/Rn2end9f/Dcx9S19lj1t+eM/f3PD27h8bn1tweN/e3L
9a2v81vsljjd6mT7nu3WCM/Y/5DRtcLzUGN/s4TlZXfWWSo8Y/+ot8v0UQ3d
OmeYsb/YX62a3lvnNHmUatSrzjlv5nbh2QnkLar07Cl3tb4XnucF8hxV7Njw
N6IjtgnPkh+ZvAieLwfyKNXDs3Zth8spwnNyIO9Sqz8ptVs/1/DsC+Rp6tnZ
x6rr+jA8XzR1rVpW/S3D2D+05uvTdd0JnjcY+0vGfhmj8wH4xf7sYz+e0/6J
3Xz+WW/DrvgOS4Xn3439809O7l03eoHwjP1HO31b5/rs2cLzN8b+2h/eOXxv
ndOcC1SVEatPaD+B588Cebh6oEDP7BeeewfydlW0+qoJyTuShefWgTxfNeh3
sYXeR+AZe6NNXgTPgXzylEo+NOFE3/zzKTybc5NKb3dwe2bieuF5AfUuq/7m
w/4hGxfquiU8DzX2J1x5brGOn/CM/cPOpoW8lrxQ+MX+k0YnjN18bp28eXbS
jDnCM/b3nfHMsS7VpgnPRYz9o/d2yrm3zmnOueqvlyYkrvKuEJ4D8eS8Ohu6
oa/OW+C5WuAcqspuzppaPmml8ByIMzlqWLTW3650/TlwzlW1ntT64ZXCM/aW
NXmR5HWBc7Rq0u1KB21PcN3+ghpg1d/aGftHxPWK1jg8NzL2VzjUt6P2N3jG
/gZlujZI6PCd8Iz9Cf3fK1shaYbwi/3pceecmTHxwjOfj43T+t6JwvNGY/96
U5eDZ1OfUfsOd2j91+wvhWf0dI3+W+nz7fnXgefDgfqJalW97evf5X8vPG8y
usSeQ9tXr5hvJzwvDNRn1NqMb67p+4LnyYF6jtra/dRBzQM8G7tUZ5MXwXN1
qb9dDqq/FTX2Nx+iE7so4Rn7f7h6NPdA5xHCM/ZPvDBqxOj4kcIz9h/+j9b9
jhKeFxn7/3Z6QvSo/PXwi/0VDhRNTc+/vsWn0yG5efYnay4K/48YvHb/0fUa
Zbj5ifl0/BN6hi7Icush8LDW4NZ69UeNwjO23h5sX1/FDNtc73p2Q1tvppqO
G3RrXL1nI626q/OfiZ+vzctz59/5/9FJsg6c9eiowNGnofsBR0+FTgUc/Y/o
RQyOXgUdgNTjjb6CvjU4egD6rOD0r+kLgtNvpY8FTn+Qvgs4/Sz6BOD0X6hr
g9MvgE+eJ3Sk8AmO7hQ+wdH1wSc4OjT4BEc3BZ/g6HzgExxdCnyCo6OAT+l7
mL4/fILTp4ZPcPqq8AlOHxA+welbSR/F4PRZ8GNw+gLwTHxCfwvP4Oh14Rkc
nSQ8y3WMrg+eZb3RocGzrDe6KXgGR+cDz4IbXQo8g6OjgGdZb/r+8AxOnxqe
BTd9VXgGpw8Iz+D0reAZnD4LPBPv0S3Ds+Q1RucsukyDozuFZ3B0kvAMjq4P
nsHRocEzOLopeAZH5wPP4OhS4BkcHYXM0Rmcvj88g9Onhmdw+qrwDE4fEJ7B
6VvBM/snem94BkcfDs/g6HjhGRzdKTyDo5OEZ6krGl0fPIOjQ4NncHRT8AyO
zgeewdGlwDM4Ogp4BqfvD8/g9KnhGZy+KjyD0weEZ/IRdPLwDI6uHp7B0UXD
Mzg6XngGR3cKz+DoJEUPZ/DSQfXDkBBwdGjwDI5uCp7B0fnAMzi6FHgWPYXR
UcAzOH1/6eManD41PIPTV4Vn8jvmCILn3xOS6APCP+uZUwieo09IYk6B34X1
6NKD5+gTkjKsOXrpqxt9NTjro635etajEw6eu09ICrPm7lmP3hWc9fY8PuvR
bYKz3p7TZ32w/pD/qIO68/usR0cHzvoW1lw/69GDgbN+vDXvL/wbXVPwewAS
kuz3AIiO0uhzRK9j1te03g/AenQmwe8NSEiy3xvAevQSwe8TSEiy3ycg/SjT
9w9+z4DmM9D3x285fzDPgn+CM/+CH4IzB4G/gaPbx6/A0ZnjP+DookVPbHB0
vHgCOLpTfndwdJKiLzQ4uj5+R3B0aPxe4Oim+F3A0fnAv/SXjC4FnsHRUcAz
5znmgOAZnLkheAZnrgSewZmDgGfpaxndPjyDozOHZ3B00fAMjo4XnsHRncIz
ODpJeAZH1wfP4OjQ4Bkc3ZTofgyOzgeeBTe6FHjmfMz8FDyDM28Fz+DM6cAz
OHMl8AzOHAQ8g6Pbh2dwdObwDI4uWnS6BkfHC8/g6E7hGRydJDyDo+uDZ3B0
aPAMjm5K9EMGR+cDz9QbmDuDZ3Dm1OAZnLkneAZnTgeeRSdu5krgGZw5CHgG
R7cPz+DozOEZHF00PIOj44VncHSn8AyOThKewdH1wTM4OjR4Bkc3JXOJpn7D
vB48gzPfB8/gzJHBs+gKzdwTPIMzpwPP4MyVwDM4cxDwDI5uH57B0ZnDMzi6
aHgGR8cLz+DoTuFZ+pZGJwnP4Oj64BkcHRo8S7/CzDnCMzhzkfAMzlwePIMz
RwbP4Mw9wTM4czrwLP1tM1cCz+DMQcCzXMfo9uFZ1hudOTyDo4uGZ3B0vKJb
NTi6U3gGRycJz2KP0fXBM/U85kPhE5x5RvgUfYqZv4NP0WOaeTH4BGe+CT7B
mceBT3DmR+ATnHkH+ARHnw+f4OjJ4RMc/TN8gqPXhU9w9KXwCU79Ez2f1GtN
XQ79GTg6NPRS4HeNbgp9Dzg6H/Qo4OhSnjT5jegzTX+fPrTUD00/mr4pOP1T
+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
/we9HWYI
   "], {{
     {EdgeForm[None], GraphicsGroup3DBox[{Polygon3DBox[CompressedData["
1:eJwtmXfgT9Ubx++55yqbkJGVVDIy27QnWjYlSUilIi1piRJpqaRBhewRoiVN
TSuRVREZZcsoJb/X+/e+f5zv87rne+/93HvuOc/zfp5TrcsdLW9PkyR5lT8R
OwN7FHZVSJLVWZLU47gUtg//nEN/Y3gBrTL9GfY7WiVdi/2KdjR8mHO/gcvB
Afs1rax+BLuEVh0uiF1MOxY+Evs97Xi4EHYRrSp8BHY37Sb4JOxd3Hc2tzkL
3k5rT//x2C208+Gq2IW0KnAB7FZaC/g47N1c+y7XNoH/oF1NfzXsGlpDuDT2
Z9op8NHYnbTO8InYlbST4ZLYpbQT4MLYIgWSZDh8AVwIfgE+D76d33qb3zoT
Lkn/SPovgYvCL8EXwqXhN+CmcDF4BHwRPJBrv+LaS+Ey9L9JfzN4Fa2uvg22
BP2vwRfDR8GjUp9/L9e+z7XnwD/SatNfHLuCVgcugS3F+a/Dl8HLaDXgoth7
uPY9rj0b/oF2Iv1FsPP1e/Ah/r8H7gHXxO6gXQ+fgN1MOxeugt1EOweujL0G
+yj37Mdc6gQ/Dj8Mt4YfhO+G28APwffAHeD+8P3wtfAA+AG4IzwQfhB+Bv4a
/hiuwTPdB3+OfY7+b+FP6W8LPwzfCz+dejznwTfCQ+GBcCv4AfgueBvP2o7j
6tjrsI/R/xD9XeAn4QFwV/gp+DH4BngI/Ch8M/wcPATuBj8NPw7345k+gs/n
ng/A8xLPkwfhjxPPgQHwl4nnxnraWVxfAXsH/TPoH8zxbxw3wVbE/kRrBJfB
rqWdBpfD/kI7FS6L3Ug7G66EXUc7HS6P3UBrDB+DvZX2D884gd+5BT4Ij4fv
p81NPIcfhj+FX9b3gT9JPD8fhecnnnvv0Ibx/3Oxz2O/07fQmHD8N3Yc5/4O
X8b/jsXex/EHib/V4/IRied/S477wX24ZrbeDx7G/ysyVydwfLXGLV+nc+h/
OF+n78EP5WPyLnwvXAGeCT8AN4Bn673g4+B34NHY1dx/Kb81nP4rOV5If1+4
GjwLHpE/02L4Hq1VeAZ8BxzhyXBvrXl4KnyX1iE8HX4iehyay+dovdH/Nn13
woXhafAj+fd6H+4FZ/AU+HY4hSfBt8EBngi/BF8FL4I/l7+A/5E/09rjPXZq
TsEF6D+gbwUfAf8Ffwgf4Jw/4aG0xTzbVfTNox2mf7/O0byg/yXNWfr/o38f
PA3ux//Wy5/p22s+wo9Fr0H5kE9pR9L/t3yP/CbX7oafpC3inCvpGwwvgC+H
Z9Ne4fzNmgvwfs7fA0/XOqF/AzwTHgJvhIfQFnLtFfQ9BS9JPB9m0Z7knE30
9c/XwgfwM/DFWvPws/Al8hfwPrgXXAe7V2sMrq1nZ449m9pnHpAfhOti39V6
4rde4NpDHD/C8SnYf7U24EYa4wL2R/KZBQrY1yi+HKY9Bp+mdaB5D9fXmqM9
BDfEluf8t1K/VwV4XOqxOgYen/ob/SWfDtfDlqN/bOoxnKMYxbM9z7MVpv/F
1H7moNYJ3AD7J02BvZbmF+c8lTp2R3ho6hiawk+kjlNDaRfAX3DPWdhS3P9Z
fevcb3wED8x91Fz4KfhCeD78InwFvAB+Ab5cvgAeBjeDv4Gfhi+Cv9Q8kU/i
/i/Cj+Y+6kP48dzXzZNfhMvA/0mU8JyPw6drTcCD4DPgL2jF4X855x1safkW
eHDmePSp1mPmePQJPCj3jR9r7snXcv5wuCz3HJN6zWoMbqX/RublLo67c1wj
f9ee9HcN1j5tNCeDtU9rze1gzdJO6ytYg7TV/A/WRFfKZ8Ll4RZaO8H66GrF
+uBYfSN8bbAmaq95HqyVOsKtgmN4V7gjXD8/p0WwjzkLbgSX1NjADYL9TSO4
drAv7AC3DPZVZ8IN4S28f0l4Jn3b4aXwM9H+spniVLCfawjXgotoHcB14JPg
G+D2wXrherht/gyXKY7ANeGWcDO4oOYpfEawzzgPPj1YX1wKnxvs2y6WDgn+
tpdIz8AnyL8n/qa6TwO4ZrB+6Qy3y59HY9s02LedDtcL1p5XKc4Ga8zr4NbB
Oq4T3AZO8nGrHzzvrlC8C/ap58CnBvvF0+C6wbGiqWJ6sK89Fz4tOFY0hk8J
9m314BODfXB9uEbwXJYemBO8LqUBPgpel4qtc4P9h/TD7GCd+yz8RLAmUFx+
P9dL0h6Dg9eEdMW7wWta8fo9uG9qLTEx2A9Jl04I9mHSe9PhIczv8zj+LFpH
309/73yO9YA7B/sV6f+pwbpvhOZJsF+U5nw7WJ92hzsF+znp7Wlw/9Rxf0aw
P5sOjw/W1y/LzwTryjfgEZqTGl+4fHB8OxYuh22ReXykQaTfB8OPBMecGnBV
7F/Bc+OkfK1pvp0JHwNfCDcO1ggXwGcF6/pu8HXBuc5N8PXBvnM8PCo4JkyD
xwXnIq8rVuba701ppFwTviL/E5xzSB/2h5um1rR9g2POBPj14Fg0EX4D7g1P
gt/UN8r8rKliKX0bab8lHhPlZf/Sfg/2z/8lvla/Oyw4Pk+Gx8AtuX51sA/v
mz/nKxz3gV9TXIPvgkfKD8N3wq8qpsAD4CnwWPg2+EX46WCt9yE8C27F/dcE
+9v76R9N/6scd6V/N/ZN+rvBe+DRcHfFJHhM9HzQ+47m+Eb6d2HfoL8t/Etw
zJXuGMM5rwXrgrHwSLg15/wUrA+lIz6j/8Ngnap18UGwvtR8mxmsm+Q3esGX
wr3gHsHx6m749uB4dRd8W7AevAe+I1hH3A53D9YRd8A3Ba8VrYsbgmPabXA3
eDnPVovjYtjm2D7096T/W44r5t9UumYm/ZPp/5XjM+V3Uq9rrdNJwbpmFjwF
7gm/AD8F3wo/Dw8Njqt3wrfmc6w3fAt8CzwMflJzkvtvw77OWL2cOedZAveA
99H/FtwZ3gqPgtfKL3NtAezP4f9hNwnYggWs76XztwfroiODv9URcNS76Dy4
EHY9rRBcJPerJ2GqBc8Lxdziwdq2Jnxc8DzlMPkn8bflZ5LDiZ9dOifDbqIV
g0tj/wiO9UnwHJREOKT1wIUDudGp2I30F6XvqOD3U06RYjfQCsPFsL/RisAl
gnV3Lbh6sMavDR8frM3rwCfA12S+Xjq5Y/4bL8PXZr7XCPi6zM/6Ctwm8zhK
p3XIPC7SG+0zj5f0z/WZ32ck3Cnz2n4tek7cwO92yH2a/GGF4LzqOPgY7H7l
GJx3cuZ87nj6K2H/pZ0AVw7O7arDFbF7aZXgMsFrsSJcKnguVIbLBq9FabOi
wWtXeUFJ+b3M14/l2d7InCf/CD+fz/Vv4VGZ8/Zl8A7OPZFrq2C3Y3fLh9G6
5HNoF+f8xfEu2g+0OfSvo39HdF1kFX1LgmsVK+HF8KDoHF/1kE9S6xblPl/C
5RLnHbqf7ruZ9nf+u7/SPk7tP+VHy+bvO5m+Ranji/KjJanHWeO9GK6aODeR
Nlyj5wzW4D/By4LrLivgRcF1HeVWX8Gf8XxFOT4YXddRzvJ17qO+h78JzmvX
wSv1ztg/aOs15/Mx2ZS4tqQc6svgnGAtvCK4zqR86ovgvOEX+MfgfOJneHlw
zelHeCE8l+c5pO8fXTdaTv8Cjt/MXJdZEZ1DKM+aH/z7m/Nn+DNxHNqaj6XG
9Pe870Di+LQ+j6fKhbek1n7SgP1zPalc79c8Fku7rUqtP5Xbrkmtl5RTr0yt
weQMVqfWtMqFH4keN41f72ifr9zwi9RzXr7rj9R6Vbl233yeKA7cqZw5cZ64
IdcAysF7RccL+d1lqTWt6iEdomOT4uo10fF0rdZ4dExch22X2Vdqrko/67up
Xqr5q3lYMB/Dvfnc1vfdY3f1/zHcl8/Jnfn4lqC9lrlupRrUOXrnxLXNc6P1
g7Re5+gYoXh+Q3Qc+RXbJTpeKE+/JTpGK0//ifOaJK4t6Bn25/Ps5uj4qDh8
o3KnxHn3/NR+Q/Fibz4Pt+XX7ciff3P+jfR9N6XW2KoNbkytjZWLrUutwfSt
f06tgVWLWJHHAn3ftan1sOoGP6bWb6pF/Cb/nbi+oTV0EdwkuIakMbkPrhcd
7xT3+mTOSbTWnssc41Vr6pc59yiWx+QunH+N/DzntUpcF12e2kep1vELXDdx
/WQ7fHbi3HNHam2voLQtdS6g+vbO1PmL6tv7UucpBfI4r/s3l4/Nn1+18V35
/Jcm3Zs6x1H9fH/q3Ee1btUhVBP7RD4cvjxxnfxAvqZUJ/8ndW6iOvl/qfMa
1WkPpc4xVZdWLfMZeBD3OZw6x1SNV4FS+ZHypH35XJRvPj1a86tOXlNzKHEt
XXNaPvMg7el8DU5NPZc1HzbIF0Tr01ex27Qm5Sui9xOWwd8pXnNcJT9Hfkjf
92TsJM4pDh8dfF/d/9vg++r+Wldvp9YGE/nzXua4otixNV8/+u3puS+dlvqd
5J+25P5pZ76+emaOi/KBAzLHV2mv8tFzQ3XXLDqvV53/yOj8WjlQgehcXrX3
ovIFifcLCkfn5qrHxug6gGrsxaNzTO0jHBFdE9AeSsnofFN12lLR+bLqt8dE
5/6qzVaKzvFV760SnZtUynWi5vzNiv/ReYrq7WfI3yXeyygTnQur3ntsdC6m
Orzypikaa3Rc2WhdofdvDb8F/4I9KlrHqqZaOlr3qj7ZlvZS4jpSwWjdqzEp
BN+cuL5aIro2ohrs0dE6WbW46tFaWrWXitF6VbWjatHaW3WYytE6WfWlctH6
XLWgqtGaXLWgItE1GdV+P0qtYaS7tsLNE9dmL4rOU6RJFkoTJda350fr4R+w
F0Zr6eXYC6K1tLTKxbThieNg0+h8Svlcrei6tGpiqgONgocz/r+nzv1V730r
c91fNaLvUutG5TgLUutG6aVLonNV7SU1i86zVmIvpT2ReA9lT+o6iWrXqsNs
yddau+gcWXXLKZlzNum99vlzKk5dF52T1svX7Pp8vXSNzmEb5jpxqsaE/8/I
nONJJ3ePzklVz+wGj0tcP7wpOj9VzbNHdExUvOgYPQ6qkV4fnVcqHrWJrgmo
pvp5as0pTXtQcTRx7f2y6LxYY7s7tS+SPwz5nJRf+jt1nUp1+HmpdZq09NzU
elUaW+OyJV/v8jPSKZrf8kXSdJWY552i6wmqu/aMzt9VP6yQxw6tp6ui6yfK
45U3aO/jHcbn5Og9KdWQp2bObxVnX8m8d/U9/G3qfEG5be3ovSTVQutE722p
znlr9JirDtwwej9LOvmU6D0m1UuXpq6taa9kUuY8XHqjUfSel+qoV0T7c+2B
jsy8r6b5/H0eZ6VXP02tXaVvx2ber1qtNR1db1G9dHTmPS3NvbrRe2SqjU/O
XBeQtvksdX6hHKRB9P6d6revZt6rWwqfl4+bNMk3qfM45fv1o7WK4vLXeaxX
Xj8xcx1BfmZM5n24VXCT6D01jee4zLnrGrhldO7/E7Z5dL1C518eXa/Qe51G
G5R4T+Ss6H067df8kLrmqf2sCZnz55+1PjLvI2rtt4quJ6i/cfSeoPLis6P3
7/SOV0bX3LTX2SK6vqFnuzqf85ozp0bvaWr/aHzm/FzP/D+cl+Wm
         "]], 
        Polygon3DBox[CompressedData["
1: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
Q//OlSMf/5Tjes/J8WnnX4l4M29d37/BS7teUqznr3Pk19ly6KHM+F7w4vBr
cIH7w0vkuH9ztIJ34/0f0+/RivN3hT+i/yPnz+v6qPnkpxLxr3uN36vwOTub
G2kf+l7Gr2H8GF6/BV6f9hB8L+Mztwl9YXxZJ0c8eqDquUGMn9SKHLhQCk1v
buxqHKI/nTbRawJ/7fygXey9yjGXndOuZWuCTeBtmqgVOuXQZuYYc42atk+K
mlKtqwZeL4XmVRurEddNobHVjmqEv+h3bkI7vJVibbnGPldr0qbRH2+e5f87
C/4Avt44B58Jvw9f57yD7696aiCfM7IVmqNJsQbUIqcx/g79q+nP3iY0/jYp
NJRa6nTG36V/jXGtTWiQf5qY42qTUxl/uwlt1pXxU+C3nL+uM/hk+E3nq+sC
no32iFo3xfnmrJxiTZnLrCmsla2ZrTUUT9ZS5lhz7RYl9MIA+t9z7tol8ttO
OfL1hiX0Sd8c+W+rEvlhvxz5b8sSemHfHPll+xLrSRNAvWO9eTq8TY58tFaJ
fLQj/DXcq4Q+2Bj+At64hD7qlyMfnQi/Bl+So54YCb8OX5qj3jgJfsP1DHdu
E2t28RQ1h2vZmLAYvGwTscKaZVl4lSZqGWPmEik0tbH09xS5wJyg1rEGWy5F
jWZtZo26QooazNq1fQ6t3rsJLWONsyq8ThO1jxp9tRQ1r9rdGmz5FDWNtdnK
JfLv1jn0gxqrlSJmqr2MiXOmiEnGyn/ga1No7kG891f4ghT3cDe4PdfgjhSa
bpi1OuNLwivAr/L+R3ltCP2Dee0MxrrBs5eIAye2ou2XogZ3zBr8qBQ1i7W5
79FneaTEZ83F992QQwNPqevLnJbrenNOOjefaWJum0NvzLH+za2b5FhLrilr
FWOKudMcaqxRI6mVzBnzw6vmiFPd6fcw9tBfnfH1GP+4jhnLjIGe45zaHz6k
ibnmHD8lxRx07htTN0sRo421xug1UmjIn+ocHZJijjp3zWFzpci55rYh8KZw
H/jHHBqlbYoaQu3yeYpYa8xVe1pjzJQiR1l7WGOoTbrV882hc6TIaeZWNf8C
KTwBawFrigVT1PTWGuY4tYqaxdynJps/RY2uVjNH9nRON5E7zaE94J5N5FZr
oo4pcqq1kjVFlxQ1h7WGNVqnFJrP2k2N3L7+XrXzJPj2FDlwYIkcaDx6uonc
qAZfKoVGVZsf3Yp8s6dhs4RmnSdFTlTL3ghfA99oDiuhIWZLUVM9UHOiXoEa
2Fx5LXw+fCm8eYma+uQUHoK1tppabW3Ong/esYQeGJyj3joKfgo+J0e9qUen
V6eH17WEx9W/zscZ6/xzPrlGXCt6dnp3emKzlPD89P705GYt4cnpzenxdS6h
EW9IUbOoHfWY9Jqs8fUmrdH1ro6v62GZHLWE81/tNLyuD2uCmUt4qHqpeqB6
rXoqelsnNOG1uL6G1vnr/F+ohGfqGj6jNmOBnp9jemh6aXoqeqd6Ypem0AR6
ZXosei16gnq3epSHpfAU9S7NiQPrenN9qgkG1fU9ZwlPwtpkWF3feop6i3qI
eqMnwwek0BRzl/DsLk6hgfXy9PzGwVc24QXunKMWM2eZu/TY9Nr0NPQy9eT0
5vQ09Hb1pa29zTHmGmsMvbXRTdQeeiB6P3qCeiN6LnpfenB6MXo0emd6Yno3
ejxnpfDU9H70mMbA45rwnvTQ9NL0oPRK9dj02syJ5kY9Lr0uPTy970Ny1FrW
XDvB++eoxcy55l49yPNSaHy9yaE5vBM9lL4lPMdLUtQkepHWDBPhSU3UEtYE
V6fwvK0V9Az1DvWo9Lb1fPR+RjbhvetZnZuiRtDLska4Cp7YRO1wJO1J+mfn
8GdOYfzAFB77vOY2eHAKT30B13MJ/+HIHP6Ev9nfPiRHvblHCb/oqBz+gznN
XHdEjlxnzjO3jay5cL8S/sxxOerTvUv4P8Nz6PsBJfyjY3Po+31L+EsjctQD
XmOv9cE5/JgDSvglJ+Xwsw4q4fecksO/Ohp+Fj4vh19yYAn/5WQTHnxECT/0
rBz+02El/NAzc/hhg0vo6VNz+F2L5qjN+vGatfQQ2lS6o3L4TSPgF+GxOfT3
MfBz8Pk5/Jzh8PPwBTn8n2PhF+ALc/hTagBz77k5tMGp/k8p9lR6Oj95fWqK
7/C73AO5DZ7cxN6I+Uv93KXmCz0YvVU9Cr0ZPZ5eVU/q/Vizqq8713ynnupR
85n6696a77rXfGWNcWMKj8baw5qke81P6i9rlm41fz5f9ZZ6Tn2nvrI+UM+p
z9S3elDqW/Wd9YI1vXpO/ac+0yOaDD/chHd0Wo7axnvivTk+h9fnnHJujfYz
UtwT7401r7XvS014UXpK96XwpPSa/vNsaj7Wy9Ej0CtQT+p9WeNZ6+kh6OVM
qfnb/PxG1ec9a/5Wn+oh9avx3bh9Bu09+tcy1r1NzCm9Vue8c+3mElpIzan2
vLVEbWRNZG10S4nayJrA2kDNpDacMYeWuqmEVlXDqmW9PhNrvvJ6OWcn1Hji
9XKN6P26pg+omlBtuDBtkRJrQq/RNe1acQ2Nr/HKtTWhxL3zHqrNJ8GbNVHT
WdtNLKGV1eBq8etLaHU1vlp/dInayxrKWmpMidrHmszaTM27EPhXCi2sRvP6
zFG1220ltJs1p7Xn7SVqS2s2a7cb4TWa0Pxqf9e03rIxxbV+Q4laQE2vtr8D
3qWJGtVaNdXPX6xqH9eoXq0xxbXrGr2nzj/XrjFFL9YYZKwxBo6t8d/YaIy8
qOa/vjX+j6v50fhujLi3zu8RdX1PrvPT+boj728Hzw7fzdhe8NwprtHjOWpM
57v6zdrTGK2XbI1n7LZmfTxFDLKWVc/rjRnjtmlFzftYihhpLXxuiVrGmtTa
1Jis92KNaKxWg42q+Vhtpga7ta4/tZk5Re/fmG2uMebpZRtzjYXGyEk1Xg2p
69X4pR4dWtfzVTXfGc/MEVfW/OX6tkbVm7fGtXY1R+gFqREH1Hw8puY/8505
44qan80l5ozRNd+bS6x3rN3dX7E2P6/Etbams7YbVSKWGlP1Jt4qkcvNseba
N0t8l9+plni7hLY155v73ylxr73naiE1r1p3n6qF3y2R+9VMaqdpJbT0fzkf
/qDEtfgvh8Pvw9c2ofHUervW7zu8avH3Smg/NYXa4pMSc8k5ZSz91DXfRAw0
Fg6s33dirQ2mlNACajq13eMltJ6aUe34WAntqCZUGz5cQouqGdQOD5XQDuYs
c9eDJXKXGkOt8XoJraYGU4tNL3HvnQPmso1aoU8eTOEPfMz4HU3kHHPPhyVi
mzFOrfRRiVykRlerX15ibVjzWfu5Rryf36RYO5eWqM2t2azdtmjF2ngiRXze
vBVrwTXj2tmqFfn5qRT+0matyM+uGdfOhSVqLWsya7ONW1GfPJTC79mkFbnI
NepavaBE7WZNaW25ZSv0wJMp/KrzS+Rac6615qatWNvWaM7PrVuR+9UA+l/G
hHbmqhSxYnyJ/KBHolcyrkStb01mbXZRidrMmtHa8ZISud6cb615WYlaWQ9D
L0OPoGMTXqDewRUlcp850Fr2yhK5zxxo7XpdiVxvztfLubpE7WtNbG08tkSt
ao1qrXpNCe/RGtla+aoStbM1sbXxtfBSTXiUepVe04dqftdbuynHHv+0JrxB
PVS91P7e4xK/2d8+Hf6ixJ6D/s86TfiXF5cY+4T2Zf089YKe+TfONfjjJjx/
/Uc9EL2Qz2hfl/CUfbbDPQy9Zu/J1Ko3/D3Wg+Ymc5T1nh65z1q4p6F3rmfv
XnTvJrz8c3jtmyaeWfD79dx9VkU/0nrS+sL6x/rGekpPz3rKesh6Q8/Oteka
tV7R43MvyXpL70+P0Pry+CbqL/d4fNbDPQT3ftzTca/YPQT3etzzcO9Dj9u9
IPfE3BvTY5+nxB6JeyX6we4lOYedy64ZvVjrPetj61vrNfdA9I/1i9XP7gm5
N69f7F6Rz+z47I57TB1L7AnpJ+uHulfknrZaSE2kP+szQz475B5VhxL6XG1g
DjeX+0yNz9a4p9WpxBp2LbuG9I7dU/JZAfdAf6kazVjgmlS7uQej3t+zCd9N
T9NnJawZ9DqNUcYq17zeps/k+GxOP9dciT0F9/rdE3WvQQ/UZy+Mcf95o63Y
23EPwr0C9yx8dsg9M/cyfAbIZ4Hco+tSwnN079J6XC/SPQT3EvY25pXwYN0b
tF61XtYDdi9Uz1Jv2Hraetb61nrXPSKfzXCP170j/X21pJpSP///uqK9fg==

         "]], Polygon3DBox[CompressedData["
1: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
A3eyV9kudpBNYf/NJP9DHo/1yLiM+A0hPNI8
         "]]}]}, {}, {}, {}, {}}, {
     {GrayLevel[0], Line3DBox[CompressedData["
1:eJwt0jkvBVEYgOGxr9cS0RINUVCKglIUVCIKOqKgIBEFpULodPwCWiL0Ejr7
cu1rRMTW2NfwnETxznMy92TmO5lb0t7b1JMQRdGwXhXWj7rUbGIU5XGBxVxn
JU9Zywc28pttzE6Kopi6rcs5yGqOsp4TbOEUOznPfi7xntd817aewz0tG2iH
K4xzlbtc416Yifvc4AE3uaXD8BwecYfHjPOEuzzlHs+4z3NOO8cNP3Shl3BG
jZlti0NcZB/n2MFJNnOcdRxhFQdYxlJ1WWewlZ/e0cA71vCYFeFcLOKtfVd8
U66ewl4VmjGHBYwxn9lhD7MYYyazmBHex3SmMkVp1klMDd843OOvkvUTftOX
jBHNuNz+/x/+AMQfSu8=
       "]]}, {
      Line3DBox[{690, 1000, 473, 689, 1106, 912, 691, 1107, 913, 692, 1108, 
       914, 693, 1109, 915, 694, 1110, 916, 695, 1111, 1006, 1206, 696, 1112, 
       917, 697, 1113, 918, 698, 1114, 919, 699, 1115, 920, 700, 1116, 921, 
       701, 1117, 922, 702, 1308, 1001, 923, 1002}], 
      Line3DBox[{704, 1007, 1207, 703, 488, 705, 1118, 924, 706, 1119, 925, 
       707, 1120, 926, 708, 1121, 927, 709, 1122, 1008, 1208, 710, 1009, 1209,
        711, 1123, 928, 712, 1124, 929, 713, 1125, 930, 714, 1126, 931, 715, 
       1127, 932, 716, 1128, 933, 717}], 
      Line3DBox[{719, 1010, 1210, 718, 1011, 1211, 720, 504, 721, 1129, 934, 
       722, 1130, 935, 723, 1131, 936, 724, 1132, 1012, 1212, 725, 1013, 1213,
        726, 1014, 1214, 727, 512, 728, 1133, 937, 729, 1134, 938, 730, 1135, 
       939, 731, 1136, 940, 732}], 
      Line3DBox[{734, 1015, 1215, 733, 1016, 1216, 735, 1017, 1217, 736, 520, 
       737, 1137, 941, 738, 1138, 942, 739, 1139, 1018, 1218, 740, 1019, 1219,
        741, 1020, 1220, 742, 1021, 1221, 743, 528, 744, 1140, 943, 745, 1141,
        944, 746, 1142, 945, 747}], 
      Line3DBox[{749, 1022, 1222, 748, 1023, 1223, 750, 1024, 1224, 751, 1025,
        1225, 752, 536, 753, 1143, 946, 754, 1144, 1026, 1226, 755, 1027, 
       1227, 756, 1028, 1228, 757, 1029, 1229, 758, 1030, 1230, 759, 1031, 
       1231, 760, 1145, 947, 761, 1146, 948, 762}], 
      Line3DBox[{764, 1032, 1232, 763, 1033, 1233, 765, 1034, 1234, 766, 1035,
        1235, 767, 1036, 1236, 768, 552, 769, 1147, 1037, 1237, 770, 1038, 
       1238, 771, 1039, 1239, 772, 1040, 1240, 773, 1041, 1241, 774, 1042, 
       1242, 775, 560, 776, 1148, 949, 777}], 
      Line3DBox[{781, 1149, 950, 779, 1150, 951, 783, 1151, 952, 785, 1152, 
       953, 787, 1153, 954, 789, 1154, 955, 791, 1155, 569, 793, 1156, 956, 
       795, 1157, 957, 797, 1158, 958, 799, 1159, 959, 801, 1160, 960, 803, 
       1161, 961, 805, 1163, 963, 807}], 
      Line3DBox[{806, 962, 1162, 804, 1255, 1054, 802, 1254, 1053, 800, 1253, 
       1052, 798, 1252, 1051, 796, 1251, 1050, 794, 1250, 1049, 792, 1249, 
       568, 790, 1248, 1048, 788, 1247, 1047, 786, 1246, 1046, 784, 1245, 
       1045, 782, 1244, 1044, 778, 1243, 1043, 780}], 
      Line3DBox[{809, 1055, 1256, 808, 1164, 964, 810, 1165, 965, 811, 1166, 
       966, 812, 1167, 967, 813, 1168, 968, 814, 1169, 1056, 1257, 815, 585, 
       816, 1170, 969, 817, 1171, 970, 818, 1172, 971, 819, 1173, 972, 820, 
       1174, 973, 821, 1175, 974, 822}], 
      Line3DBox[{824, 1057, 1258, 823, 1058, 1259, 825, 1176, 975, 826, 1177, 
       976, 827, 1178, 977, 828, 1179, 978, 829, 1180, 1059, 1260, 830, 1060, 
       1261, 831, 601, 832, 1181, 979, 833, 1182, 980, 834, 1183, 981, 835, 
       1184, 982, 836, 1185, 983, 837}], 
      Line3DBox[{839, 1061, 1262, 838, 1062, 1263, 840, 609, 841, 1186, 984, 
       842, 1187, 985, 843, 1188, 986, 844, 1189, 1063, 1264, 845, 1064, 1265,
        846, 1065, 1266, 847, 617, 848, 1190, 987, 849, 1191, 988, 850, 1192, 
       989, 851, 1193, 990, 852}], 
      Line3DBox[{854, 1066, 1267, 853, 1067, 1268, 855, 1068, 1269, 856, 625, 
       857, 1194, 991, 858, 1195, 992, 859, 1196, 1069, 1270, 860, 1070, 1271,
        861, 1071, 1272, 862, 1072, 1273, 863, 633, 864, 1197, 993, 865, 1198,
        994, 866, 1199, 995, 867}], 
      Line3DBox[{869, 1073, 1274, 868, 1074, 1275, 870, 1075, 1276, 871, 1076,
        1277, 872, 641, 873, 1200, 996, 874, 1201, 1077, 1278, 875, 1078, 
       1279, 876, 1079, 1280, 877, 1080, 1281, 878, 1081, 1282, 879, 649, 880,
        1202, 997, 881, 1203, 998, 882}], 
      Line3DBox[{884, 1082, 1283, 883, 1083, 1284, 885, 1084, 1285, 886, 1085,
        1286, 887, 1086, 1287, 888, 657, 889, 1204, 1087, 1288, 890, 1088, 
       1289, 891, 1089, 1290, 892, 1090, 1291, 893, 1091, 1292, 894, 1092, 
       1293, 895, 665, 896, 1205, 999, 897}], 
      Line3DBox[{911, 1005, 685, 910, 1307, 1104, 909, 1306, 1103, 908, 1305, 
       1102, 907, 1304, 1101, 906, 1303, 1100, 905, 1302, 1099, 904, 1301, 
       1300, 1098, 903, 1299, 1097, 902, 1298, 1096, 901, 1297, 1095, 900, 
       1296, 1094, 899, 1295, 1093, 898, 1105, 1294, 1003, 1004}]}, {
      Line3DBox[{251, 474, 1106, 252, 488, 280, 1211, 503, 295, 1216, 518, 
       310, 1223, 533, 325, 1233, 548, 340, 1244, 563, 1150, 355, 578, 1164, 
       370, 1259, 593, 385, 1263, 608, 400, 1268, 623, 415, 1275, 638, 430, 
       1284, 653, 445, 1295, 668, 460}], 
      Line3DBox[{253, 475, 1107, 254, 489, 1118, 281, 504, 296, 1217, 519, 
       311, 1224, 534, 326, 1234, 549, 341, 1245, 564, 1151, 356, 579, 1165, 
       371, 594, 1176, 386, 609, 401, 1269, 624, 416, 1276, 639, 431, 1285, 
       654, 446, 1296, 669, 461}], 
      Line3DBox[{255, 476, 1108, 256, 490, 1119, 282, 505, 1129, 297, 520, 
       312, 1225, 535, 327, 1235, 550, 342, 1246, 565, 1152, 357, 580, 1166, 
       372, 595, 1177, 387, 610, 1186, 402, 625, 417, 1277, 640, 432, 1286, 
       655, 447, 1297, 670, 462}], 
      Line3DBox[{257, 477, 1109, 258, 491, 1120, 283, 506, 1130, 298, 521, 
       1137, 313, 536, 328, 1236, 551, 343, 1247, 566, 1153, 358, 581, 1167, 
       373, 596, 1178, 388, 611, 1187, 403, 626, 1194, 418, 641, 433, 1287, 
       656, 448, 1298, 671, 463}], 
      Line3DBox[{259, 478, 1110, 260, 492, 1121, 284, 507, 1131, 299, 522, 
       1138, 314, 537, 1143, 329, 552, 344, 1248, 567, 1154, 359, 582, 1168, 
       374, 597, 1179, 389, 612, 1188, 404, 627, 1195, 419, 642, 1200, 434, 
       657, 449, 1299, 672, 464}], 
      Line3DBox[{261, 479, 1111, 263, 493, 1122, 285, 508, 1132, 300, 523, 
       1139, 315, 538, 1144, 330, 553, 1147, 345, 568, 1155, 360, 583, 1169, 
       375, 598, 1180, 390, 613, 1189, 405, 628, 1196, 420, 643, 1201, 435, 
       658, 1204, 450, 1300, 673, 465}], 
      Line3DBox[{265, 481, 1112, 266, 1209, 495, 287, 1213, 510, 302, 1219, 
       525, 317, 1227, 540, 332, 1238, 555, 347, 1250, 570, 1156, 362, 585, 
       377, 1261, 600, 392, 1265, 615, 407, 1271, 630, 422, 1279, 645, 437, 
       1289, 660, 452, 1302, 675, 467}], 
      Line3DBox[{267, 482, 1113, 268, 496, 1123, 288, 1214, 511, 303, 1220, 
       526, 318, 1228, 541, 333, 1239, 556, 348, 1251, 571, 1157, 363, 586, 
       1170, 378, 601, 393, 1266, 616, 408, 1272, 631, 423, 1280, 646, 438, 
       1290, 661, 453, 1303, 676, 468}], 
      Line3DBox[{269, 483, 1114, 270, 497, 1124, 289, 512, 304, 1221, 527, 
       319, 1229, 542, 334, 1240, 557, 349, 1252, 572, 1158, 364, 587, 1171, 
       379, 602, 1181, 394, 617, 409, 1273, 632, 424, 1281, 647, 439, 1291, 
       662, 454, 1304, 677, 469}], 
      Line3DBox[{271, 484, 1115, 272, 498, 1125, 290, 513, 1133, 305, 528, 
       320, 1230, 543, 335, 1241, 558, 350, 1253, 573, 1159, 365, 588, 1172, 
       380, 603, 1182, 395, 618, 1190, 410, 633, 425, 1282, 648, 440, 1292, 
       663, 455, 1305, 678, 470}], 
      Line3DBox[{273, 485, 1116, 274, 499, 1126, 291, 514, 1134, 306, 529, 
       1140, 321, 1231, 544, 336, 1242, 559, 351, 1254, 574, 1160, 366, 589, 
       1173, 381, 604, 1183, 396, 619, 1191, 411, 634, 1197, 426, 649, 441, 
       1293, 664, 456, 1306, 679, 471}], 
      Line3DBox[{275, 486, 1117, 276, 500, 1127, 292, 515, 1135, 307, 530, 
       1141, 322, 545, 1145, 337, 560, 352, 1255, 575, 1161, 367, 590, 1174, 
       382, 605, 1184, 397, 620, 1192, 412, 635, 1198, 427, 650, 1202, 442, 
       665, 457, 1307, 680, 472}], 
      Line3DBox[{277, 682, 1308, 683, 278, 501, 1128, 293, 516, 1136, 308, 
       531, 1142, 323, 546, 1146, 338, 561, 1148, 353, 576, 1162, 1163, 368, 
       591, 1175, 383, 606, 1185, 398, 621, 1193, 413, 636, 1199, 428, 651, 
       1203, 443, 666, 1205, 458, 685, 686, 687}], 
      Line3DBox[{459, 667, 684, 1294, 444, 652, 1283, 429, 637, 1274, 414, 
       622, 1267, 399, 607, 1262, 384, 592, 1258, 369, 577, 1256, 354, 1149, 
       562, 1243, 339, 547, 1232, 324, 532, 1222, 309, 517, 1215, 294, 502, 
       1210, 279, 487, 1207, 250, 473, 681, 688}], 
      Line3DBox[{466, 674, 1301, 451, 659, 1288, 436, 644, 1278, 421, 629, 
       1270, 406, 614, 1264, 391, 599, 1260, 376, 584, 1257, 361, 569, 1249, 
       346, 554, 1237, 331, 539, 1226, 316, 524, 1218, 301, 509, 1212, 286, 
       494, 1208, 264, 480, 1206, 262}]}, {}, {}}},
   VertexNormals->CompressedData["
1:eJx0nHdcze8f9yMKlZ2SbEWlgZaoF6lkZWRERBlRpEFGopA9U0ZTUUo77b33
Os1zTpGV0VlmIuM+3Xfvfg/fx80/1+OTnI73ua7rPV7PV1NtDq7bPVBERER+
sIiIqHAVc5xtmT2FDVrNh504nj2lE7dKnGa0O7Bx3oUx6txlJo4OnGzU7tAJ
g+YFZinpbGQnjZhu/q4ZX/MNHVLSO7HztJxYjUgrGva3r3fXaYLCfVfFGhEO
TqhO3xW/pBXXRqSJnzzWAF3n2rPxSzjQHjZO+p1nK9JVA6M+xDKw+Jmk9jtP
Do5XTJs7MLMVe3q05sxMqMVRSSvzgZkc5Ba9XHCD34p3NZbbbvAr+1YONrl0
KCvLt8FxjdGKGbJl0J3uu1JZnovrh1d6OS5pA8/nWLeKRhHMtFUeOS7hYsin
tD3Dd7Wh4+7tX7Pn50GyPOfC8F1cjDZx35Pp3oZHy32Diu3T8OtFyIVMdy60
nzxbqn+jDV9sDgTtXxeP9hPj7fRvcKFeZGG9z78NKq4HT9Va3Mfqj/fP7PPn
Yure1Y8Vgttg4qFxzVzR2sBgRECVQjC3P44UV4rnIcMh/jHNzP71o+qVoprc
Trx+/z0qfRILVbvXtklptCAsVCpgtxQHXu0qb3V2srDC+UvSRrcmxKtLhMza
xIGclsfdpAcs+Pm4yQ1Ib8C8dzxd/wAOPskZWma1szDr/sqGyPcMnC1+MGnx
Uw6uTF344ek4Nuo+7jk2engdMkTbt96X5SJ86KX50aZslMt0Oy4xq+pbuRD4
7jfLO8TG/G9HNhfbl8F1/+PzvfGJCVh9O9OPDdGkF4FiR4qwWfpmfmg4Fy5W
BS3LMtgwWaRX9sA5D3X7IuyCyriI9AxXnN7ARsdIk/KOgjT8CfYYbv2Ki4Lv
PgFSHWwoZPvN+HolHg/GnVkX3cWF2wmf0mcCNuYyV1UE3roPT51nHusG8rBO
pCT7xGc23i/6EyS9y9rALmhPgo4Yr39f0j6l/UlxpLhSPLNYvJXth1r6V6mY
fcvjtnHgGjvb3SKlBZLtzRVuWU0oMvwzcPRjDuacDekY8qUF3a7+B8W7GmBS
MsN92CcOxpX+kPRQFX6OR9VHT5lWj/S5YR7SWlxwxRI3rbZh4n1tg/GlxXVY
czZtYKIzF/NmXvoj4c3EKNP9TjsfV/WtXJTx5ymszWLCaWb50rVZZagyvjYx
qZWL7tGqnPznTPjIXN+/MrMIDbss43r/v7oRZuck/jChV3HkqfjjPFxrU7h+
VYmHYV5stbfjWJCbOu7u23HpmKSsp7TTmAeTFJe2CbNYKJh/adve4nikGjV/
PryZh6fZJ9JOzWXBvebPzI1V99GmbzjZdzcPFzU6NL9rseCfwvZ/FGxt8OC4
dGvoPl7/OadzT+ed9iXtU9qfFEeKK8VTdomnyOkfTf1rzbWIhwO/c1B74OKA
X5rNqDbhN6loNOLmg/nhveeX62G20sO+Ga6Ps0O+bKrHr0eb+YwLXKScUu66
H9gME1nlsAfOdeAsitzXu9923OTd/VYh/P66hpPrBlb3rTzsjJ+/ddrnZig5
e71YOLQccQ4Twp20eHi2iD+0VFp43hzOXzn3swgxk78Pv2jNw61SKTnfuS1Y
7vUD3c/yUDXhtd7cczyYm+bcX7qsBZ15usY6tum4dlTSWfoBD69lzg55tLkF
noNXRNr0xON35JwuTioPlt0NZy7tbMETrWErfQaG4J7NH6ZnEQ/uV8Kms/a0
oGDYdvn7TGsDxQP3rFLLeP33Jt2jdH/SOadzT+ed9iXtU9qfFEeKK8VTenv6
AnHbxv71XYXi3t77bcDGRLW5fo0I2CZdWnKpHqGsFEF7CxfeTmbXTcobceu4
vcngkDoo/v68pHw8Dy4HJju8+NyI3RGuLx02VPetPLRqLRYLG9+EUUFim76v
LYfbzd1nV13iYXfsjiI9vSasq/n1Vk+vGPKNsl9747N3xaClGzY04a3iDY1R
o/LRarWy/elTHgrnswwT7Zrgf6ohaeGTdOj4KLiU9PCw6cJsLfnjTXjj89DO
XS0Bm8+5lmwZxcemL8G19mea8FB99tMPWiEoXMZn7JzER8inyGEXzzWhJKrn
lPdwG4Nldbs/5k7j9+chykuUj+jepHuU7k8653Tu6bzTvqR9SvuT4khxpXjG
LJ+nr5VR37/uqG3QipjKw+yEGkZ7Rz0e/1HX2FBVh7UB0yb+3/M4T332PqkG
1Ai+Oko/qO5beXhuGJhwVKMBMYNFtnT5luOJ15htwSwe+PbXj2SZNaArS+UH
z6UYe5eHHXoyhA+lY1UaHNsGBAe8jeAhH/qi6xLa1fn4EFAgUuLWgKjZ81K3
/UiHeKz28gIzPgruhE6Wu9QAxXuOV9dYJQBj9iyW3C18ndNBegXeDei+FPt7
4J4QcKOY588483Fjo/TJR74NWJxse3nkYhuD90eP+A86wu/P65TnKb9THqK8
RPmI7k26R+n+pHNO557OO+1L2qe0PymOFFeK5+vJjov1vtT1r6tybXt2R/FQ
99qodZssA/bD4kwfdlb3rTy0cSQiDugysPTQlGsPG8vxR0Fe0Lt/dst8cJy4
noEJXgfezYgshl5X+rV55nzMGvpFMMueARdvt7FK9vmw+CQZc95N+PUZwZc2
uDNwpvl+lMjCDLx5l5+5wZ+PhZIxz10uMXD7wpDVhy4loHiydf3yRD4g0SJn
7c2A76VKeRvvEHwsVutWyuNjxmCZd+K+DCgbVq+MPmhjoDB5R3t1Eb+/Tvq7
buL05/W/8zy3Pw/9nZe4/ffm3/cor/+c/33uef378u99yuuP499x5cH2x++n
7io1/1n5KLo/f5Hk7hrs+31o8VGpir5nPpjKgW5Ft2tgKfjCtmUV9z3z8STs
W5dSXg1EJf4YWHvn9z3zUTee8YjxvAZ/zJNdrx/J6HvmY+S4rV27umtQEfMV
ZxIT+p75sJ4bqfJQrBZ+D9Z/fpQR0vcswLi4etO9ErXoFmzX/3HPxuD/PQv6
606qQ6n+pDqJ6iaqlyivU56n/E55iPIS5SO6N+kepfuTzjmdezrvtC9pn9L+
/DuONf3xFIuP38xFRf8atn+0+PtgPkz3O2zzsq/AkNTms0aDSnBhqph2ST0f
231WuzVer4DGly473eR8HH3tdmbaLz6OVc8+cDy6AscCClKPR2fg1TfjsmuT
BLh7bcBY5fwKiL6TdpJvToCBVvuwBG0BDM8ufvy4ugIqKZF5zs9C4Lfh2Iil
xgJwI0Y7sBgV6NTIuHkr18ZgSNqiZbOWC/rreKrrqZ6nupPqUKo/qU6iuonq
JcrrlOcpv1MeorxE+YjuTbpH6f6kc07nns773/uypn9/UhwprhTPa3cbitco
lfSvaQ/2KffutzvhvitUTUuQzws6uakhH0ayZzqYcwSolygZe2FHCZY+rZyq
z86AiFL43Y/mAlQ/7kyZ5FwCk1Fn0ju+COO839L05n4BwgpqTPNOlOCmpJy3
7a8Q7H756biLmwDDdLyazTxKMDFjQpjscxsD+0/mix08BP19EfVJ1B9RHU91
PdXzVHdSHUr1J9VJVDdRvUR5nfI85XfKQ5SXKB/RvUn3KN2ff5/zmv7zTvuS
9intT4ojxZXiuXrakjmO3Pz+9ZVnw8SibQIoplpNgUgBBhiZNVuJZMLHdfCi
3viovvJKvDq8ADFLDX9fHZ4I1YGBuUfuCeNsVxI7T6YAXnZ7vTVlQvEo45Ou
X4QAcsyXTrJyBRCRFrOd12Nj8G275bXHMYL+PpP6Tuo3qS+iPon6I6rjqa6n
ep7qTqpDqf6kOonqJqqXKK9Tnqf8TnmI8hLlo7/vzZr++5POOZ17Ou+0L2mf
0v6kOFJcKZ5WxXf9j07K7F/DlfavufJAgOjtViXamplorgu2rp2WiFOromd3
Zgnv4SHycZaLMzFnZTduqYRC/E3MG1QJYGw88Ua+cSayS9fF3Ru50+DUp0kH
rBmC/r6d+njq36nPpL6T+k3qi6hPov6I6niq66mep7qT6lCqP6lOorqJ6iXK
65TnKb//nYdq+vMR3Zt0j9L9Seeczj2d9/592bdPaX9SHCmuFM87eGTgNCex
f5V70e7g1SzAprL8pmV6iWhyiHRYrheK49c+uoq8EeD2+23e0Bd+HpdNLEdM
2WmglfvKQ40j6J+D0FyE5iHUt1MfT/079ZnUd1K/SX0R9UnUH1EdT3U91fNU
d1IdSvUn1UlUN1G99Hder+3P75SHKC9RPqJ7k+5Ruj/pnNO5p/NO+5L2Ke1P
iiPFleIZv/z3d/Ulof3rH+U3St5fBQg6ITFC2SQUql0F4/Yq7TTo2uRcU/pD
0D9XojkTzZdoDkJzEZqHUN9OfTz179RnUt9J/Sb1RdQnUX9EdTzV9VTPU91J
dSjVn3/XSbWgeonyOuV5yu+UhygvUT6ie5PuUbo/6ZzTuafzTvuS9intT4oj
xZXi+eNPbuqF2TsNaO360/tHALfjfmpBJ1n964UPzWYWjp1QZTmc7NjEwqav
g7dqCOtn0eazbhqmnZB8eFXuQUUrZMecfjHu6A0YFIddmFHPxXWt3ORrla0I
dj4/+5CxtUFnY/ZT2QZu//fRv6PvX9Lzwtxuyw0EjGk8pB3UhpEBk+9oBwn7
mrudLRFLvfvXkVOm+O78JTy/HuF2eSo7DVgbTEtsTb0R6DEhJeG3oP/90ful
90mvSz+HXp/eH71fep/0uvRz6PXb1Wa93V3H6l87Rr+Vvyp8fflEEfNYrTYE
jK0tSdvngXxlnfE/07j9z/T39HVWbHhUuenp/lXKpKjJQhj/3VsZNQVZLMyP
n3/ymiYbQ4eFBF7T7MSPK7Hs98L37zVb6mf9Fg/YcX/8fC98/yoa579/0rI2
KF63JvG98PVffHl+RTedi4sTXx/4LnzduneTQzcKP9+lj/z2RQtfn16Xfg69
Pv17ej16Hfp59PPp59Lr0s+h1795mbHzUjUTyqsNFbuGsrDuXFbAt+ROWC6z
qz3xiI3FtosCrqxtwWgrS8GJR51YPFuyVVWKhcT1beMV/VrwbC3DyONHJ9pl
bpWWvWHjWYuK8zp2E36I1b0te9OJoc2nkhWWsxBTJtL2cHIzvCQ1M2Yt5EBO
Z9LRrmmtuGF7XFZ8cyO+LXlwt2saB4+v2jkqnmXh5ORS1Qy/RmxXnBkVfJyD
g8ZRg9ItW2GqJH1zUE09zDYdUki35MAuYsmZ1xksXNb/dDBPogHGIUm5t5I5
mP94n3vRjVaI37UakTqfgW2rnIOKbnDgLXtLRZrPwrus4tkKLgywj8w7nMTl
4FHuqfSa/FYk6w100ZtfC+eZheya/P/N89//Z56f1jRj6cGJbGxJeFFfH1qL
42+XDH0y+X9z6Yr/zKUXqh5f2zSuDXvVhu40ulUGr2PTXZrGcXFH5Zrg7kE2
LjHP99i9K8Pmleui1I5yMfwlw8FYrw0MwYi34crFmBgedsNYj4v76U8ePb3J
htfAp7pph4rxYXGq9rhALsIOJ+2OsBDu1/SIhrUp+ah5EH0+woKLa1Gjwg3j
2DjdJXuz4Us+Qu4lS+VkceFcMr6qyqkNMlIHH3w7mIVaaQ1+lRMX7Xt+M2+X
siG7P1R3SFMWXltZWnxs5mKT86X8C15t8NvDSnhvk4xm27RXF7y4SHGTD41n
syEyuer3prhkxL5NNzrD4WKaweSAjFttws9DMOaidzQ2alvlZNziIle84ZzD
OzYk9TYOGfU9Gs7mH4tW/+DCRn+m47KANgTP8wv87R+AQUccvZcFcNHZcqMw
9QMbSwIe7M8fGIjDAaWqC0V5eJpT+XP9zBa0b7a4uW0dC6NfLLyjK8oBd/Sy
HfK7WjD1wveo8vYWKNxapD11OQevfvkV7whqgRYiBqjbNOPItvToCdc4uCIl
VsxsasEMmfVhjaxGGHd0bayv4YB1elGN3jAmpjQYq0gvbUDAlG5rgQQXQwUa
H/UXMDGgXHXp9scMPBwenDzRmAtDi/Rkv71MlH4LUdMeVAeOlMuHocf/Nz8f
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
Uv8PC37Yvg==
    "]],
  Axes->True,
  BoxRatios->{1, 1, 0.4},
  Method->{"RotationControl" -> "Globe"},
  PlotRange->{{-1, 1}, {-1, 1}, {0., 0.9999998571428621}},
  PlotRangePadding->{
    Scaled[0.02], 
    Scaled[0.02], 
    Scaled[0.02]}]], "Output", "PluginEmbeddedContent"],

Cell[BoxData[
 Graphics3DBox[GraphicsComplex3DBox[CompressedData["
1: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
xj5XLthygY1vOMRDTHzDYf1u4hsO603wDBzrIxPfcPDnwTNw+J/gGTj8JfAM
HPM7eAaO+Qg8A4f9NPENft/BM3CMTxPfYD5NfMP9/rUzCjbX32PjG25I5uc/
KtvSk3nG3xFfMvE6h3iIidc5rN/BM3CsN028zmF9BJ6Bw58Hz8Dhf4Jn8Zu8
vwSegWN+B8/AMR+ZeB3bT/AMHO+7idfx+DTxOssn499Xb33n7G92nbnouEe7
VZdsD3A9q+6Gjpcu2M5xPeAndL5gYFH1Dke/4xoUHnPkG1+9+Sa3C3x4dt3Y
Mx4e3Nu078z9jN+0639KurWOT9O+q/t88Pz6+iv4d3zOHjvLHlq4pSu3j/Zs
P4Gb9vk5op+m/9z+50MbBrS2F+A65NHfn7B/7Q5Hv+MaHHf27GEvfvmiM/c7
g/P1xydmTt9Wfxr3D/935iMfv1BStJPvw9+H7f+vI48vLOb20d+dlX94Y9Nq
6T9/7l4HXL9ky6G9zf3cvrk/MO1w/83n8v+hfT/6mwN9LSjY0DqYZf3VHLzQ
pbA4tPOzQ5jXX8DL3EvRtY7xVPT/y93G6Fpt7q9wr/gr4x9Gt1W513x7jK/3
/0/2cGop8O9bwp+N7hn6fOAv+n645+hzgEe/FmfcWtVec+B/bXAv0P8Bvyi8
vcsWN5P+rvuzzaXoCnxAeFvruKPfGbe/a16FZ1zPOq3qqSWt8w541niCeQbe
99IFy49v9bvAs74/yTwD/314e+u6AzwDvyFqZx3zDNx/rxrmGXjvqP008wy8
9TOvCeNO4Bl41HyveuYZ+Jqwl9c0Ms+6P1uYZ+D1HUpfD/MO4BX4QM1/oMep
8Ay8MvrcD5hnjQfMM/Da6IOXMs/6/lXMM3Dfz0rmGXg6+sNa5hk4xjl4Bu75
2cg8A791ykEjwzgeeAbeEH3uZuYZ+E0Ph/9Rzzzr/jQyz8DBP3g245F51u+9
8Ay8+6TwE6bw/2u8lHnG9S9RPxczz/r+Fcwz8Oj21nUxeAbuv+8a5hm4f6/X
Mc/AT4za38A8Az+gX0hEmnnWn7uJeQbePHj/fmFeADzr/tQzz8CraPyDZ+AX
k/0Bz9qOCs/An4xevPeYZ43P5XaBZ6JxspB51vd/wjwDj26/PMk8A2+gcQ6e
gc8iew6egUdft/N65hn4tZ07bQvjzOAZ+KHR56aZZ+DgHzzr/mxmnoHD/oBn
8x4zz3peEp6Zj+h9mcQ8a3wW8wzcv6cL+PP0/UuZZ1yviXhYyTwDvy1qp5x5
Bt41mtdWM8/AN5GdAc/A32hlM4zbg2fg1xH/4Bk4xj941v3ZxDwDX0r2HzwD
70bzL3jW87zwDPw37U9qfSPfYp41PpN5Br6a7DZ41vcv4n4Af4rGOX7HtYrs
OXgGXkw8gmfgvaL2VzPPwME/eAY+nMY/eAZ+Etkf8Kz7s5F5Bo75FzwDH2X8
GHwPyzNw7y9NMH4d8OnGr5ta6uepecavw/1lxq+bWupvW2L8OvgXy41fh2vS
+HVTS8G/9uumlmL8a79uainsj/brppZeTPZf+3XoT7Xx66aWDiT/BzwD3yO6
Mv+B7q/wDNzP4+OYZ41/wDwDT5M9Ac/6/oB5Bt5A9hw8A68jO4N+az9uOfOM
61oa/+AZ+O1kf8Az8Eay/+AZ+N9o/gXPuj9rmWfgf1b+ZzPjteq5Ngf6+QvP
wP17/RrzrPH3mGfg3k+bwTzr++cyz8D9fJRgnoE3Ev/gWftrS/j7AB9O9ge/
49qX7D94Bn4Yzb/gGfgj5P+AZ92fSuYZ+B/I/wfPwOfQ+wKe9fskPAP3du8l
5lnjbzHPwDHOwbO+fybzDHw48Q+egado/INn4KPsOCb8DLL/+J7A8V7gd1yf
pucCnoGfTP4neNb9WcU8A3e0/gLP2n+Udam2T8IzcN/+88yzxscxz8DryD8B
z/r+D5hn4Fka/+AZ+GayP+AZeC3Zf/AMfB09F/AM/A7yf/D99ecu5N9x7U/+
P3jW/VnKPAP/C61/wTPwtJpfmgNt74Vn4GdF42QY86zxl5hnXqcS/+BZ3/8W
88zvJdkf8Ay8lp4LeAY+muZf8Azc0fsCnoE30vMCz8Cfp/cIvAA/ndZf+F33
J2CegZ9H8QfwDPxVNe8h7rMtxjPw30X9v4d51vgw5hn4Vhr/4Fnf/xLzDHwM
2X/wDLyJnhd4Bj6O/B/wDPwieo7gGfjn9H7x+w+/jJ4vv2+En0frX/Cl+/MB
/47rxRT/Ac/Ax9N4AM/aH4nF6xg38TrGTbyOcROvY9zE6xg38Tqjk+B4HeMm
Xse4idcxbuJ1jJt4HeMmXse4ideZfjbHcIonIa4U9I/s2JvMP3Dcj/gfcB8H
kHggcO+X7uDf8X9+PF+RJ54o7eP6v9Sf3HFH6T/ae4HGubmf4p/Sn3z84P/G
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
mo2rY/+Ujatjvw94Bm73xwGHbt/G1aH//6F6deirbVwdemAbV4d+FTzbPAV4
tnFgG1fHfl4bV8f+RxtXx35JG1cH/zaujvFv4+rYPwWegWO/j42rY3+Kjatj
P4WNq0P//0P16tBX27g69MDg2eYpwLONA9u4Ovbn2rg69vPauDr4t3F1jH8b
V8f+PhtXx3408KzXwfG4Ovb7xPTqNM5tXB37KWxcHfr/H6pXh74aPNs8BXi2
cWAbV8d+cxtXB/82ro79vDauDvtj4+rYL2nj6np/n8TVsR/NxtWxf8rG1bHf
x8bVsT/FxtWxn8LG1aH//6F6dcsz/o6rjatj/76Nq2O/uY2rY3+0jatjP6+N
q2P/qY2r6/2SElfH/j4bV8d+NBtXx/4pG1d/ydgX4NifYuPqr5rnCtzo/831
h+vVrQ4cuNGfOxuHt3F1o28PbBz7h+rVoSf//xpXh64ev9s4/w/Vq0OXbuPq
0LHjd1zz6dUR3y7R/cobVzc6c/N/sfvzxtUNzv2vpvZNf1gP739v4jrbL5Ts
/2pRNcZfluvPfHTw4NqWFsFx//ZTz9rc0oI8ZobvX/Pm2vTuOO6/+INBJxZ3
yRKe5vubonYaGcf9k3454f4uhXhPUnz/nosm1Le01DOO+/d6Z/WikiLoLpJ8
/5FX7pttadnMOO7/S5d9Di4sxvuZ4PtPaSja2tKyiXHcP398r77VJbDPCa7z
Mz/ip5px3O+/F/KkSae/71rGcf/E6PtC15fi+ztF/a9kHPffGPUf9ijN9/t6
qqsYx/2nRfVUV+K5mP4sZVzzvwTP3bQfMK7bn49xZe7/gHFdv/19HoeoBzZk
+75XVJds4fHG9b0jHhp5XAGfX7rPQzJO0ox7fjI8ToDv+ezebxVV1/F4AD7x
i/a/Ke4izx24fy9q+fnq/shz1J+7mp+Xbr+cnwvwqqie7XLmX7eziHnW9djn
MZ8an8x8oo6ab6eR+QR+0X2VTfIeZRg/5ZNxH8t7lGbc19GtYz6Bt7vzN090
KdzEfAL3/dnIfAL342oD86n7s4b51J9bwXwC97zJOAfu6wB/zHzqdsqYT33/
HOZT4+8yn6g/d1Hf9itKiuqZT+D3jnnjM7EzGcbfXnv64WJn0oxXHlp+poy3
FON7RnZsI/MJ3D9H5o3x0sgurWM+dX8qmU/9ueXMp25/BfOp66UvYT6BV0Tt
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
8+BNzKfnuYrO9Uszn7jOvnre46MvrzV8rgvclX8/vNO2ajM+1wd/brr905Ii
G3+rpLgor0OJh9XB4qNrjzytahXz6XleE8x7Zejdlzd/wnx6nquCTO9+P5qU
Xsx8ep7XBud2umvLiZMSzOecEVVjXmzt57uX3XtQ1VMSfzvlhp+sPLa1n4Pe
+2z87vG3q8pntLT608G8xvJdJUUSL7rwnuP7D6xeEjz/Rp+Lqp5KM5+/Pnzc
itZ1S3DT+e9d0f6kjcyn91OXB2fXV/Ua0quG+fT8rwxmd+58dvuTqplPXD/5
YugtUw6S+JvnudzGkYi3RcHJUZwwyXx6fpYEvznknYsevnUF8zk74uHj4MKe
jy2s77CM+Zx4Yp+CktZ+PvnM849J/YeMe7l8Qf/1rf38uuuUlonp+cznsHtO
X/l9y6qgzz8mDZ3WW+Jvc6P2y4Of7Dqmw+7xt7N9YCA4/6zxN4u9ygb9Llj1
/S8KJwQPPH7Sa1MOqmE+f17+9qtziiYHF/95/c7jh1QznwMiv2p6MKNqfP8F
y9cxnxdFz6UkKGkc+v3g/XkeIf7nBVcUPzar92kSf8P12bC7hRJ/u/Pw93sc
29rPeZH9kHiRxycEU/Ya/L9LvpF40XH/fGLlzNZ+JisfPmRIr8XMp79/enB8
VeKLGy5dwHxeesOrYQ43uOPwe7s3D57LfP4qamdesCOb+vLZL6cxnwduGfja
dy1BMOOcK57aPf5GuOsw+fUa8TeyAbXjqhYv/qeMn0zwe/+57pQH9h9x6YK1
zCf103Vbdvuo/fuxfxLQ93IfftBQKvrDJN0/wW08+ZF9ynkeSRBe7CZccsrd
1SUcR+LrYdGfJV6E86j8eXmLmU96jq7v4AdvEnuYdvTc3SvnZVYO3v8j5vNn
fpy42n3PvrZX+1nMJ40rd9Rt4bmTkl87149D916nA75saRnHfNK4dfVHjLlS
4o3ZgMa5O3X/698XPyQTvOTfC9fx5TO63b3/auaT3iM3Zs5picubOS8T0Hvn
Zt538YPit7S+5/49dbevevDqmx5OMp/0XrvBj7/4bpdCib+RHaD4mMSLcC1u
bL/qms5lzCfZGbfi6M/792ofMJ9kl9z69RfOPH7IHOazp7dj7k/f7PdkqBsH
n2T33M0F981+asl7zCfZSbdjWphIl/gb2VVX8uTLK2U8ZAOyw+65Sz8vmXJQ
OfNJdtsVTL19wW5+CPFT5Tpd9d068YdTAc0L7pZHr797SK+VzCfNI25r9eeX
bjtvOfNJ844bd0digMSXEN9YT3E/iRfhvEx/buB85hPXTf22793v0jnMJ82b
rnN0DulM5pPmWbd4W+X05cdPZT5pXnbvjV73Ru/TJP5G87jrUZhcsHv8jeb9
Vntw7IiB1Z8wn+QnuNVdf/L4tby/IEM8pF3y5Z5zhvTiPAvx3Do/9Xqxxwr2
Q1LEc627+aAuzZc3L2E+yc9xLUfs2CHvKcbnBrfiyp0nSdwy4ciPovjPHMPn
Jvf5Sz0uk/k35chPc7vGZifKeOP8lbtpvwNb/Q2OMxDPtW7w3N8Gxw95l/kk
v9HtPPHMv4vOp8mRn+nuT5/6zu7xN/JL3eRF//NgV86PZInnejfiX2e9dC2/
Rxniud6NW/7yzgXLOY5NPGRc45FdR/ZjO4/1+GZ307htA05iPwTxzM3uuCkV
P7uN12tYH9W5zx49s0F0AhifDe7HD3ccJfG0pKN1gZsanff3AfNJ6wh3+PSl
Iy9dMJl5pHWHu/Sc8JzQScwnrtf8Y9URpa9PYD5pXeNufvastbIfpMnROsjd
v/OEn+0ef6N1k3vp7F9WzS5iu0E8b3VPDFx96Qz2EzLE8xZX3a7pksNKZzKf
tI5zr0XnvnFcMaB1n1vQr3/loXw/4huN7k/dPr5hOrePdXqju+SK1wbM4v5g
fbTVPX3As3268nyadLRudQ99WLLvLh4/KeJ5izu7+0crR18+nvmkdbGbsvH8
0w8tfZ15pHW0e7Bgaeb1oaOZT1y3z2331pJvXmA+aZ3u9vh5eL6hxN+euuz7
PqE+6cCnf/bXvSc8znzSet992vLh6pXcTiag+IDLzplWN4E/N008b3dv/OiW
dp24nyniebt74LutY8bw90oSz03uP2c/+tuvmAfEi7a5iy5rrhfesH7f6ZKN
x43+Jce7EM/c4f64/ovaNdc8ynxSvMX9a0zFr4eNfoL5pPiMe/D+8BzDYcwn
xXPcC8/WnPQE358lnpvcmrfDcxUftbwFG56+r+G+OduZZ+AP77i/4LQq8Qfw
fxee0b7bpLSs9/F3g/N1+a+qVi355k7bvhu+YN+FLfWnMa7jWj04rgX9z6RH
D9zQ0sL95P28gyKB2vv83P0VOjde1/D9Xv83n8eDvl/yCLjf6zwlj6Dvl/lX
t495NmXulzwj7vf7JiTPqO8X/1z3x/rnuF90Hbjf61Hl/Hp9v8SXNB6vf6PP
s0iY87M4bmDaYR2d09+X4w/mftbpmf3irCMy93Pc0ml+OL5k7uf4J98fdaeL
2B19P8dRje4Y7TSZ+3cwrvdhQ/fYxPuja6K/c16Mcf9/nOdl3PPGeXZzv+QH
gfvvJf6zbkf0AMD3UM9D6h9An0/3m33kHN9g3OurWY9k7mf9ktP9kfkVuH6O
Waf7mWWegftxuIV51udhsf9jcPbnGfd6WF5PMe75/Jh51vUJJK+t99OzvoJx
3/8q5ln3R+I/wD0PrDti3OvJNzLP+n6OUzn9vTiuZeorZJhn3c8G5hm418M3
Ms/Q0Xsd+0TmWdcP4Dw74/77in+o71/CPAP344TjDKYd0WMAbxf9vpp51u1L
XAi4Hz+iswXuv+8G5pnroyr+ZV+dHv8Zp/tfxzwD34PuA8+6vkU986zPw+L8
r8HFbwQ+N2ruI+ZZ3y95cOB+PhVdAXD/fHl9zbjnWfbV6fYlXsT1SCN+RO8K
/E7Ff4rxl9X4l/1zc5X9yTjd/zTzDNzzXMc86/oKm5lnXV+B/XzGR6n9FBnG
a5TdTpv7Oc4Q6DoZsn7U7ch6HPgexBd41u2Ljgu4tvNSt0+P/xTjlyj7I/vk
apT9l30wg5Q/kHW6n2nmGXihei/gN+I8rPHMs8Ylbw7cj7fZzLO+P8E8a5zj
PwZfxjzr9lcwz/p+jjs5jUt8XuOi39Y469AYb6fuyxic45YG38A8a7yWedbn
YbHOweCsM2H8F8qepM39ovcAru150rTDcQ/G94x+X8Y86/Z5/PM+qrnK/si+
+UuU/U8xPkfNv7IfS9uljKmPwvFMp/uZYp6Ba3vVxPuj/ee+wjxrfCLzDNz3
n/U/ga5/w/Efxj2fpcyzbkd0SsD9eFjEPOv+sE6M96Xdqe0/48ep+Tdl6pqI
fl73h+OWjHdT652s0/1czTxz3deo/TXMs65PU8w8a5x1PozrcZ42909nnoHP
UfwnTTvzmGfg7YhH8KzbL2Oegev5V/Zn6+eSkvqlyv+U+p36fZH6fHeqdWjW
6X4mmWfg2o9q4v3Rm6K/s56H8bSyJxnGj1P+SdrcL7om4H3U+E+adiROBXyv
6PfZzLNun58L75ssVf6P7Fvto96XFONztf/vdH9Eb6/rx8j+Jt3PZcwz8MHq
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
zDN0kp5njmeyTtLzwbo+1kna+Bt0kg9F+fppzDPyHNFj5/qQopO8XM0jopP0
52mKHhV5F697kf3gmCd03Qn7E9dJ+ucg8TfYca17FJ2kH58p5hm6E6+rFD0q
dJL+Pa1mnqGTvEDphEUn6e2u7AeHTtL/WeJvyK8cpHUprJP0umip/4/1s9Y9
ik7S778QfTW916Q3k/M+sD73vCXteObzu8Ev4kpkX5lnXM9X4010kn4fBNdH
ZZ2k3l8gOkl/v+iroZP087jUN4BO0rcncU7oJG38DTpJ7Q+ITvIE7RexTnKI
0lmJTtLrNlkvyjpJr1fhODzrJHW9CNFJWh07zZs0H0r8jeZZ92204YrPdeKr
143XMs/QSeL8a/AMnaTXAcp+Aegk/xnxI/U6oJP0f5c4J3SSNv6GeIrW8YpO
Uttb0Uka/SfrJP2+LamfSX4XnXdcbsbzBuJZ4pzkbzhbf0brJDcYnjeRjpfP
I2OdpD+3WvSRuHp7IvE36CSviOyh7H+BTtKPQ6k/A52k75eNc8bjb9BJer2r
xN+gk9wr6qjE38jPJ7st+7kQb+2g/CLRSVapeUp0kn6cSj0l8p+Nbll0kl2U
PRSdpB8nsp+L1k3Em8Tf0H+cPw6ecfV+tezngk7S8yz1lKCT9M/Hxjm3xuJv
0Elq/1Z0klF3lsf2J7oT6fxx8AydZLPyo0QneZ7yi0Qn6fmT+mC0HnQY1+CZ
1uk8TsEzdJJe7y3xN4oDkH2L7U90en+W6CS9fyX7E3H1PEt9MOgk/e9ct411
klrXJzrJJrVuFZ3kJjpPHDxDZ7hA+fmik/TziOyr9f1vov0aEuf0/G9juwue
oZOEvQHP0En6/vK5mayTbFL7DUUn6e0tn7PAOknfT9lXC52kfy7CL3SSnucY
n8ExSncqOsmFRscInWQ/rXtknaTVQ+J6m9ZVsk7S6jD9FXFCrssXEB4UqB/R
SWqdUts6Sa1za1snqXV0beskrZ4q9/2ik7xQzYNt6yS1HqZtnaTWfeHn/5xO
8lu1zmpbJ6n1aW3rJGequETbOkmjM8xzf8b0R84/yn1/1rQv5yLlvr/J3C/n
JeXW9YlOUusARSepdYPpILfOUHSSF6o4legktT5KdJJa19e2TlK/R0mXW4eZ
MvpJOecrt84z43LrQrPmfjnnK7fuVHSSWteXDXLrAEUn6dcBkhfLrTMUnaTW
JYpO0urKcuse29ZJal1l0uXWYaZcbt2m6CStfji3LjTrcutIm1xu3anoJLWu
T3SSWgcoOkmtG0wHuXWGopO8kNZ54Dm3jlF0kvp5ta2T1O9d0uge5VxF/bly
rmJunWfGtCPnKubWkTYZXOq6a92gnPsJ3Psjlcxzbt1gOsitMxSdpNUtAD9G
xVVEJ6l1j23rJLWuMuly6zBFJ+k/V84D1e3LeaC5daFZc7+cR5Bbdyo6Sa1j
yQa5dYCik9R58LTRT0od6dy6RNFJFqrnIjpJrXtsWyfZouyb6CS7aDvmcus2
RSdpdBpGPyl5tNw60iaDyzka/mp1fVYnyXo2g7P+jXH9XKxOks+XN7isTzUu
dQP058Z1klpXaXWSXJfD4FwvxbTD877B7f4yqyO1Okk5/yW3ri8b5NYBik5S
6zTSQW6dYSrIrUsUnWSJjp8EuXWPbeskNf+ik1yn5wuXW7cpOskXlf3PuNy6
0KzLrSNtcrl1p6KT9OtBOacstw4wE+TWDaaD3DrDVJBblyg6SV0PWXSSWvfY
tk7S/y71GXLrMFMut25TdJLazmRcbl1o1uXWkTa53LpT0Un6OIacr6f1inK+
Xm7doOgkjY4ryK1LTAa5dYyik9TzRds6Sf+75B1y6zBFJ6l1m6KT1ONcdJLa
notOskTNm00Gl3PioOvznyvnQubWAYpO0uiyjO5RzjHMrUsUnaTmX3SSWvfY
tk7S2H+je+RzrBjX/o/oJPW+mIxpR843zK0jbTK4nG+YW9eXDXLrAEUnqXVu
aaOfZB1XkFuXKDrJQjX+RSepdY9t6yT1/Cs6SeP/uNy6TdFJan8mY/STci5n
bh1pk8HlXM7cur6s0StuYZ5z6wbTQW6doegkq9X4F52ktj+ikxyg7H/bOskW
5f8kje5R8pW5dZuik7xDrb8yph3ZN5dbR9pk2rHny9g4ntVJNjHPGt/KPGs8
yzxrvIF51niGeda4nJOi8Xw6ScnvaDxlrla3aXWSSdMvqwu1Okmujxrgmq+e
ZL2pkwOd5N/VfnnRSb5m6jBAJ/mA2o8vOsnrlK5SdJKrov+L6yT930fEdJJb
lP8vOsnPVP0B0UmOVfUKRCf5ldaFsk5yv+G710MQneQjqn6C6CS1v912PUld
50F0kj/NWQ9nS6DrSIhO8gBVd0J0kqVq/hWdpK23AJ3kYPXeiU5SjwfRSeo6
G6KTfM7UW4BOskiNB9FJnq3Gg+gkc8c5EV/icc46SV2HpO16ks7UK/PzVn2g
66KITvLPph6O9zMaSJf1IfPs+W8w50SITlLXdRedpK4DIzpJUzeJdZLPax0d
6yTXq7o0opPUdTlEJ2njnNBJ6jo5opPUdXVEJ6nr8LRdT1LX+RGdpK4LJDpJ
8AiePf/5dZKz1bpAdJJG58Y6SV0HSXSSum6S6CT/qeosiU7yZKVzFp2kjXNC
J+k/n/VprJMs1Xo21knqulKik7R1fXHNV09S6xxEJ+n/vph59uMzv05S190S
neQBqk6X6CSHq7petp6k6KvJ/6P6tFJ3AjrJUlVnTHSSNs4JnaT250Unqeue
iU5S10kTnaSuqyY6ydNUHba260l6Pj8x47mCdVLgmda9VK9Y4m90X9Co9Ves
k9yl646yTtLoYVgnqevgiU5yk66HwzpJG+eETlLX5ROdpK7jJzpJXfdPdJLP
ab0c6yR1XUHRSeo6hG3Xk8ynk9S46CTLVB1F0UlOV3UXRSep6zSKTlLXdRSd
5N5qHhedpI1zQifpnz/vD2WdpK5LKTpJXcdSdJK3qrqXopNcqupkik5S19UU
nSSNX+ZZ//zwepJGF8Q6SV1HVHSSvU29MugknapTKjrJQ5VuX3SSNs4JnaS3
Q7w/lPNPus6q6CQfUXVZRSep67iKTlLXfRWd5CZTd53sgKkT2LZOUuvG264n
qeviik5S19EVnWQnVXdXdJLanxGdpI1zQidJ/hfzDJ2krhssOklTZ9vUk2Qd
F+skr1V1jEUnqesei07SX9canvPrJHUdZtFJ6rrNbdeT1Hpv0Ulqey46SV13
WnSS/vtJnBM6Sc+/xN8QrzQ6K9ZJ6rrZopPUdbZFJ6nrcotOUutsRSfpfxe9
E/kbsTgndJK6rrjoJHUdctFJjlR1y9uuJ6nrootOUtdRF51kiVpniU5S12kX
naSu6y46SV0HXnSSum686CT1+BedpK5LLzpJzXfbOkldJ190krquvugkdR1+
0Unquv1t15PU/rnoJHPHOaEzlPgbdJL6nALRSXZS5xqITvIadQ6C6CT1uQmi
k9T2X3SSvr2NzDOtB/PqJO25k9BJ6nMiRCepz5UQnaTW1YtO8jy1zmq7nqSN
c0IneanS14lOsp86X0N0kuXqPA7RSS5V53eITlKf9yE6yTPU+SCik9T+T9s6
Sf8eVzPP0Enaeuzo/7XqPBTRSdp6v9BJTlfnrYhOUu/3EZ3kVUYPCZ3kh0Y/
CZ1kT103knWStxg9JK756kn+29Sx9Fe7Dz2uk/Tnl+P841Y7V7L/q2HcG+dz
A8c53f5c8wbGh2zf94pw/O6IzjWvY/zGLvscHOaDbo7ar2b8or7tV4T5OH9/
BeOnFO85PMyHeryM8SM+3uN8yTvgB9+jFvfz90sWfXtf+N5R+4yPr/56Z6iH
of4wflfhl7eE/hT1n/GhAzcHYb6Svi/jx1y84YMwX0a8Mf74Je+1jrFpTvOZ
DbafetbmMF+g+cwGHx08uDZ8zprPbNDpyn2zoR3SfLbO5/dVNoV5NM1nNrh3
zBufhfOG5jMbvD3vtl1hHlnzmQ3eivCNhk/ENTcYPrPB9ee/sl/4/mg+s8Hs
4d0PCePYms9s0GFl2f+E76HmMxukj5neJ4xjaz6zwejJtX8O32fNZzYYe+th
DeG4lfPj0c9l479acmdvzXPr+v6DQSeG+RfNcyZY8+badMi/5jkTzC/d56Fw
PtQ8Z4JTPhn3cZiX0TxngrfXnn54mJfRPGeCjvWr/hbmZTTPwGsMz4hfrjc8
Z4LGTiePDfMymudMkOx546/DvIzmORMMevOof4Z+veY5E5y19rF1YV5G85wJ
juhU1nlDyXzDcya47YD2HWtK3jc8h/189YR+xw0yPKeDSb+ccH/Im+Y5HTRF
47zR8JwmO5AxPKeDjj8bMDGc/zXP6aDy0PIzw/dd85wOOkX31xqe08HIH/+6
UnjBD+KUKcNzOrj7s4+7h/615jkdNNaNuja0A5rndDDznXkHzmpdJ2me08HY
EX/oH657NM/poF3Nkkdnta53Nc/p4PULam6e3Tr/ap7Dfo4/ZGXZAMNzKtjr
ndWLwnGoeU4Fey6aUB/aAc1zK/7s3m+FfGqeU0G7O3/zRJg31zy33h+1v9Hw
3IpH928wPKeCPX7z+jniT+IH8cgqw3OrP7v0po/CvLnmORV81/n5/4TxAc1z
KvjvB1fdG647Nc+p4Gv38y/CvLnmORU0jVu+Mcyba55TwapHP5w+qtWf1DyH
/Xzn29/ecLvhORn8JRqfdYbnZHBkZIc3G56TwcQv2v9Gxmcd42v+vl8i5FPz
nKTxv8HwnAyayva6UeapMsbv/Uu7qbI+wg/ijqsNz8lgxBvf/UHmL9EJdn9g
lwt1IJrnZHDurM82hut4zXMyuOrfTY+G85rmORms/Wsy9JEMz8lgz7kjWwa2
ro80z2E/p9a/99UthudEMH98r76hndQ8J4JTGoq2hvOU5jlBfkit4TlB9meD
4TkRlEbtrzM8J8ierDU8J4Ktz9z68911ZcD9tcLwnAh+tfXEhaGuSfOcCK5b
fcOToZ+geU4Ewz96tjC0D5rnRHDrZfd1HFj9keE5EZR1XTMu1DVpnhNBz3/f
3Ofi1vW+5jkR7FEw+9OOxTdbnp23w+ssz25+5G9UW56d9zfWWp6d9zfWWJ6d
9zcqLc/O+xvllmfn/Y1Vlmf6fYXl2fWN/I1llmfn/Y1Flmfn/Y2E5dl5f2O2
5dm9EvkbUy3Pzvsb4y3PzvsbN1q74SZG47DK2g3if621G87bgdXWbjg/D1ZY
u+H8uC23dsP5efBTazecnwdt3U6s25dZu+H8PLjI2g3n58EF1m64D6N5cJ61
G87PgzOt3XB+Hpxs7Ybz8+BYazecnwevt/Mg8VZh50Eat5V2HnTePpfbedBV
RfZ5lZ0Hqf0Vdh502yL7/LGdB523z0sMz4j/LbLzoPP2eYGdB523z6V2HnTe
PpfYedBdHdnnaXYedN4+T7TzoGsX2edX7DzovH2+2vp17rQzUn1Cf9L4de6F
HRPXhe+v8euIz+WG57S7r/yMHuE8bvw6GudLDM9pV3rrPzrKOGS/zvWb+emz
wiN+RH9j/Dp3wei/NpTw+NztfPgTF68I+TR+nft08S9mhHkG49e5LnM79AjH
rfHr3JTsupYwvm38Onfrib3HdC0stn5d6zUxa3j3y+06xfn5a4ldp5DdWGrX
KTQ+FxmeM8RnmV2nuIrIbiwwPOP+j+w6hezGPMMz4jWz7TrF3RPZjZl2nUJ2
Y5pdp5DdmGzXKWQ33rLrFLIbY+06hezGCLtOIbvRx/CcpfE83667aTwHdt3t
vP82z/CcpfE8x667ic8Sw3OWxvNMw3OWxvN0wzPr9O26250fjefJhuesezoa
zxPtupvG8wS77qbxPNauu2k8F9t1N43noXbdTeP5PMNzk+sR8fy+jQsRn5MN
n03E57s2LkR8TrRxIeLzLRsXIj4nGD4R9xpv40LE51gbFyI+X7FxIeKz2MaF
iM8RNi5EfA61cSHic7Dhs4n47NHb96sR3yvodOC6vUN94/CoXxJPa3fe2rKw
nsCY6L3ZwPjW+9c8Esbvj4n6Vcl45dTKM8J6GleVz2gJdQbA5zWW7wrrmZzt
HQDGzz9r/M1h/YS5I6rGvLhb3f36I8ZcGa7vTrnhJyuP3a2uT8mTL68M849X
HN6tQOpcNAUHbDl2xMDqT9zBjy3sKLrXpmDyov95sGvr81iRvbFHGB8H/tLZ
v6yaXTTDPXXZ933COCbwA5/+2V/3nvA4PbfNeI+C9r/Z1iust3DwDc+MCXU/
wHtcV9I7rHfxWsTPasb/+tDj54X1Ri685/j+A6uXMP78G30uCuu93HjBqu9/
UTiB8QceP+m1MJ/reRBdyXOXfl4S1u969V/DG17cbX/M6q4/efzazp8QD5Ln
HPGvs166tnMZ8cB1LoNhA1dfOqP3TLc84kF0r5+2fLh65Tcv0HOvhX0OCvfp
elBY3+On0fddw/j8DvOvD+ur9Dx83IpQFwL8pvPfuyKsb3NU+duvzimazPjF
f16/M6wvtOGCbSsu4voNmaD85Z5zwrpw/ntJnnDc8pd3hnX5/PeSfEt1u6ZL
DiudSd9rB+MNc6bVTRg6mvpfhXk8aLj/ySfD+jD+vVnO+Nn1Vb3CvOHA6P2Y
zvjMqvH9w/pIvp+yP6DxyK4jw7pPvp+Sfxj73V0rN9XPpH5uZ/yNH93SrlPp
68TbSvhvwezOnc8O6wVdFI2HEsZnNQ79PqzXVBN9ruiIF/TrX3lo6/f1n7ud
8Qe/2zpmzOXj6XvNo3aSwZXFj80K9S2+HdijZPD12Y/+9qs333THRHZA1muF
jQfcGeaZvd1oZHzIfmvvCetyeDsj67K3ad3h+dnAeOXvix4J69J4O1bJuI//
VND9ixg/OYq3J92dh7/f41i2Awk3L7K3K8hOs+7GrBc4r+58fDLhTo3mATkX
0Purc+h+WX/9+OGOo8SP3cr48AOe7dOV/StZfyUbjxv9S7af4CfpSrolLg7r
mfjvu5nxQZ27bgrryVwejYcaxo874pn7wno+nofVjC8+uvbIsJ7SJ9H7uITx
3xzyzkVh3s3zM4HxKXsN/t+wnpjnR/J7B04f91WoC/K41O34/KUel4X1oPw8
Wc/41J3PXVl+zQduZ8TPVsYf+rBk312t48Tzs4PxP67/onbNNY+6oyIe8H1T
rmLKwtvCujq+/7WMn/T1sGPDukaehzWMz31l6N1hXanZkX37mPELez62MKzr
ddw/n1g5k+1GyiUrHz4krKvm52Gpo7NrbHZiWG/Qf1/Or7rDpy8dGdYj9d93
C+PndP9oZVgP1n/fHYy/P6bi18NGP+HOmTu2f6jL9HjanfNF99dCHaP/XmxP
XKZ3vx+FeuCJJ/YpKBF74p565vnHQr2cf15sT9xxVYkvbmi1G77/Ugf6D+fU
fBrq34ZE/Wd74v618fzTD221G77/bE/cQ/e/8eam+mHu3Kifa9l/PrfTXVtO
nJRwL5Uv6L+e7UzGfd11SsvE9Hz3xxteLRjAdibjBhx+b/ewrp3vj9RFfrBg
aSasr+v7s53xF56tOemJVn6G3XP6yu9bVrE/2ecfk4aGusqTo+c1j/HmbOrL
Z7+cRu3D/mTd6rc7LyxvHT/7bhn42nctAftR08+54qnQj6L5l/2BDU/f13Df
nO3Wb2n9+/ndE/+5q7e53z284/6CMP9o8obkb39g53c3fMG+C1vqT+P6F2hf
1yVoDuRz/Q/qOOj7F3DdAeCo54B98vr+ct7XDRz75bEPGbjerze1FDj2I2Of
J3DobvQ+ySbShW42+/qaSN9Rz/u+gGM/HfYp6f5s4X01wH2z2wyfWVN/o5nx
3Hxmub6B5jPLdTA0n6Kb03xmKZ9fafgEvtbwmQ2wj1vzCXyj4RPtbDJ84ntl
DJ+4v97wCbzR8AkcdqHZtN9keM6Y+j/NjOfmOcN1WjTPGa4ronnOcB0MzXOG
6zZonjOmzgB4zgTYF695ho6mxvCcIZ1X2vAMHVmd4TkTYF+n5hn9aTA8Qydl
9zWBt62G57Spj9HMeG6e01yXQ/Oc5rofmuc01xXRPKe5DobmOc16Q80z9HRr
Dc/QjVYbntOss9M8p9meaJ7TrP/SPIu+TPOcDvQ+RN7PxvOG5jll6k01M56b
5xTXEdI8p7jujeY5xXVaNM8priuieYYOscLwDN3cGsNzivcFaJ5TZIc3GJ5T
ZIc3Gp5x/ybDc4r51zyjHbt/D7w1GJ6Tpt5RM+O5eU5yPRnNc5Lrcmiek1z3
RvOc5DotmuekqSsCnqGPqzQ8J1mfq3lOBtC/aZ6TrA/VPCeZf80z7x8yPCfZ
/miewVvG8Jww9YuaGc/Nc4LrJmmeE1yXSfOcYHuieU5w3RvNc8LUaQHP0MGV
G54TfA6d5jkRoG6D5jnBfovmORFgX7zmGf2pNTyjHbtPFbzVGZ6xz36C5d/p
+m/x+nHmubjczwXtT3fmeTn4Ofp5YZ+0tf8Jh3pZ+jmifWuv0H5gnm9BwYc5
3y/sZ11onjv0jMvteKD9CEvNeMB+2VV2nDjUk9HjpKAAdVHM+OG6KHr8FBR0
yT2u6FphxhV0r1V2vNF1tRlvBQWww2YcOthtPQ7BT7Udnw7zqR6fBQXwW8y4
pfvht9tztWPjma4brH12ur4Z2+c84xP3T7b22cGeGPvsMD6Nfea6HMY+m3op
bJ95/Bj7zHVvjH3mOi3GPvP+EmOfuQ6Gsc9ct8HYZ1NnIHaOkfU3nK7fyP5G
Hp5x/0TrbzisH42/4TBvGn/DoW6Y8TdMXRT2N4jnhdbf4Poqxt/gujfG3+A6
Lcbf4Loixt/gOhjG3zB1G2zdgdXWf3Z6HyX7z3l4xv0TrP/M49z4z9gnZ/1n
hzpsxn82dTnYf3bYd2/8Z67LZPxnfi+M/8zPxfjPXKfF+M9cV8T4z6YOhq0H
kbTrQaf3BfN6MA/PGbYzZj3o4LeY9SDzb9aDDnXtzHrQ6TpsvB50qBtm1oP8
XMx6kOsy2TpKeF/MepDr3pj1INdpMetBU1fE1ulYZuMbTu9z5/hGHp6zPP5N
fIP5N/ENHv8mvuFQJ9DEN8w+C45vONRhM/ENrhtm4htc58rEN7guk4lvcB0h
E9/g52jiG6ZOi62fkrBxOYe6DSb+xvyb+BuPfxN/Y/5N/M1h37GJvzm9T5bj
b7yv08TfeB+iib/xvjkTf+N9Xib+xvuSTPyNn6+Jv5n6NrZ+zVTsw8H34roi
qH8BHHUwUK8BOOo2oL4AcN//1Q774YHD78L+beDYx419sMD9vFbnsG8TuPcH
6lHvgHHsN8S+ON3+FtT7YNx3A/FJiSdj/xHqfdA45Lofd3r/jnHUqUA9BeCo
qzDP23uDV7nxflwwjn3rFxJPwLEv1dZBxD5K7PcDjn1/2J+m8S0O+6l0O9vo
uUueF3U2UA8C+FKqC0H7chlHHYO0H4+MY9+951vyvNjviX2JwLE/UdcVywTY
TxfdvlueF/u/6LwLzOPBZNo/hXoKwFFXwbcjeV7UAfD9lDwv9kVi/xvwa2kf
H/ZrAce+M9Rr8Hgq6El1G1BfADjqDGAfIHDsB/TjXfK82L+GegQeTwZ3U10C
33/J82JfG+oZeTzhsmTH1is7kHCoz1Oi7ECC98d3UXYgQePpI67TARz1Oroo
O4B9cku5DoK+f7k9z5n3lRcrO5DgfdD6XJME79stVnYgwftMS5QdSPC+SOzf
05+L85klz5uhOhXVyj4keR96kbIPSd5n2kXZh6RD3ZISZR+SDnU2XlD2IelQ
F6JQ2Yck7/suVPYhyfuUi5R9SPK+2heVfUjyPlDsVwSu9y0iDp5yfaleSomy
Gyl3JNW3KVJ2I0X7c+dy3RDg2G+OOhfAUe9C240U77MuVHYjxfuCi5TdSPE+
1lHKbqR436V/vpLnRR2eEm1PHOrGoL4JcNQ5eUHbE4e6HMae8D5lY09cR9pX
a+wJ7wP176/keVGvZoCyMxmH+ipTlZ3J0L7vD3m/LnDUR9J2JsP7TFGfxeNZ
V0fjP6nsT5bG1RTejwp8GO1LRZ0Rj4sfhX2qhJv9quy30N938P5VtIN9rHup
+5rcqWofK+dt6XfE4ST/i/2qNv+r64O1nf/N438GefzPII//GeTxP805jJL/
bcntfwZ5/M8gj/8Z4Bwxm//N439Sf2L+J+Xp3o/lf/Osj/LwmQ3yrI+CPOuj
IM/6iPMaNv+bZ30U5FkfBXnWR0Ge9VGQZ30U5FkfBXnWR0Ge9ZGZPyX/m2e9
n4fnTJBnvR/kWe8Hedb7QZ71PvOjeUb+N7beD/Ks94nP2Hqf1gWx9T7dH1vv
c50bm//Ns96n72HX++kgT/wqD8/pIE/8KsgTvwryxK8472bzv3niV0Ge+BWf
K2Hzv+DH5n/JX4nlf1FnyOZ/88SvuI6R5hm82fhVivJcsXhsHp5xfyweG+SJ
x3L9Kpv/PTriIRaPDfLEY4M88dggTzyW7GQsHhvkiccGeeKxQZ54bJAnHkvf
w8Zjk0Ge/EIenpNBnvxCkCe/EOTJLwR58guBrp9s87+x/AKPN5v/9d2I5Rc4
v2zzv/62WH6B8542/5snv4DxEsv/fpcz/5g///tdznxiIsidNwRu84OJIHce
MMF8ap5xv83rAbf5O9mvavO//mrzcfhcm3cDbvNrwG0eDe3bfFn+/O/RNK9p
/u3P//v879E0/n9o/jd3Hj9//vdoei9+aP43d94/f/43tx4gf/7XPxerE8if
/82tH8if/82tK8if/82tN8if/82tQ8if/82tT8if//V23uZ58ZMv/2v1DEmX
R5+TZ3wmXR59joNO1eZ/8+hzOD9o87+6Hr7N/8b0OVxvzeZ/i+h52fyvHw8x
fQ6tH2P6HI6T2vyvbyemz6Gr1ecgbxvTm+XhGffH9GYuj96M8302/7ue5neb
/8V8pHlG/jemN3N59GZ8DovN/1YTnzb/6++P6c34XAab/9XnaNj8r9WbSd5W
85w//5tHP+ny6CddHv2ky6OfNOf7SP4X85HN/+bRT/I5ODb/W0jj0+Z/wafN
//r+xPSTLo9+kq5WP5lxefTAeXjG/TE9MPU7pgd2efTALo8e2OXRA7s8emCX
Rw/MdtXmf/G+2/xvNY1Pm//Nowc252jY/K/VA2ddHn17Hp6zLo++3eXRt7s8
+nbiLaZv53Grecb9MX27y6Nvp3Zi+nbCN8fyvx6P6dsJj+nb6XNj+na6Wn17
E+tqbP43z/4Ll2f/hcuz/4J0L7H9Fy7P/guXZ/8Fz0c2/+vtamz/BdmT2P4L
Gs+x/Rcuz/4Lc16Mzf/uiOV/EZez+V/EkWz+F3EPm/9F3NLmfxHvtfnfFloH
2fxvF/LPbf63hPwfm/9FfNLmf6tpfrH534FkD23+d5eKA0v+90k6b8Lmf/tH
8eo3Y/nfZ+ncBJv/9XnV92L53yrKO9j8b9Sdp0pi+V/kHWz+19ftnxfL/+I8
d5v/3Urnj9v87/VR3urdWP7X1/OfGsv/RmmB0TNi+V9/XvmsWP7X56fKYvlf
n0dYFMv/4tx5/95I/jdN56Tb/O/pVK/Y5n+Rj/PtSP73CcpTdKH5Fvg/qE64
zf/ivHub/03Q+ew2/+vxIJb/BQ82/+vrppbH8r84593mf/9G9Zx9/yX/i/Pf
bf73VVoX2/yvz4POiOV/B5I/0EXZAaz7Fsbyv4Xk33ZRdgB1fz+N5X99++Wx
/C/qVRcrO4BzIdbF8r/+uiGW//X81sbyvxiXNv+L8+5pvgLPDufmVCv7kHRN
9B4VKfuQdPq5SP63Inohl8fyv0dFiX85DwC4z0tWxvK//w1f923Vsfyvz9/V
xvK/3m6kY/lfn9+si+V/m1RdaMn//jLKOyN/Lfnf+ZGOZUks/+vrcq+M5X8/
o3Fu879et7A6lv+NzENpOpb/vSJ6Yepi+V9/LsbmWP73jei9rqfnK/nfA+i8
khJtTyj/viqW/51POh+b//XvdUxPQt83pifhuug2/+vf98ZY/tefL5CM5X+9
HVsdy//6z10by/96+9wYy/963rKx/K/PF6+J5X8viPB1sfyvz8tvjeV/R9H7
bPO/x6i6xDb/OzWW/11I+WKb/71O1SVuO//7e6X7kvzvoJzrJumPzf/WkL9h
878Xk39i879zyJ+x+d9u5P/Y/O9LZLfp8zn/u0d0HR/L/9aS32Xzv3PIT7P5
35fIr7P5303kB9r876vkN9r875sRPjiW/8U8ovmUPKbN/+K8Dpv/vZjiLTb/
C52Szf96PmfG8r/6/CzJ/+5Bz9fmf2vp+dr8L56jzf/iedn8b5qei83/vkr8
2/zveOKZeGO8gH5s/ncUrdM1z5LHtPnfrhQ/sflfPw/a+BXOL1oQy//6czk+
iuV/a+g5ap4zQbvoOjuW//0FPS+b/72E1n02/4vnYvO/x9F7ZPO/GeLf5n+b
aZxrnjN5eEa9+ZWGZ8lj2vwvdG42/wt/yeZ//f/beGyaxvOiWP53lOJR8r9+
PCdi+d9aei42/4v3yOZ/cR6Ozf+m1TpL8r/wM23+dzyNc81zOg/PKRqfFYZn
yWPa/C/8SJv/xXlANv/r218Ry/8OoriBzf/iuWieU8Tzolj+F++Lzf/6cV4a
y//OIf5t/rcPjX+b/51P9sfmf5eRndc8p/LwnKTxU2V4ljymzf/inBqb//X3
VcTyvxfTc7H5X/99P43lf2vofdE8J4nnZbH87y/ovbD530uIf5v/xfi3+d/j
yP7Y/G+GxrnN/zaTPdE8J/PwnKD+rzc8Sx7T5n/3IP5t/tfjVbH8r7erq2P5
X4/beCzaSRqecf+KWP53T+Lf5n/9/VY3AjwRy/96v212LP/b3vcjlv89MLqO
Nzwn8vBcUODtc00s/+vxtOEfPxaXfJzN/26g72nzvxRPieV/L6b1m83/+vdu
Qyz/S3GNWP4XumWb/8V5TTb/68f/2lj+dxQ9P5v/1edPSf43zzhxHrd53oKC
Wnp/bf53Lr3vNv87h+yqzf/eSXbY5n9fonZt/hfzo83/pnO+7wlXSn6Lzf++
SvOdzf/+kuyDzf+OzznfJdxp5If4323+Fz82/xuzz/Rc6sz4lDymzf9eTHEB
m//1z2VjLP+L86Zs/vdOGm82/4vx6X+X/K8fD6tj+d9LyP7b/O9xNP/a/C/G
g83/9iF7bvO/H9G8afO/H9P8qHlO5uEZep2M4VnymDb/i/OjbP7Xv1+bYvnf
GuLf5n9fJjtg87+XkN3QPKcc7L/N/86h+dfmfzeS/2Pzv33I/7T537nk/9v8
72W0/rL53+tpnat5TuXhOU38NBieJY9p879+/Gdi+d9ZxL/N/2L82/zvXLLP
Nv9bo3iR/G87+h42/wt7a/O/GOc2/wt7aPO/OLfQ5n8zZK9s/reZ/BDNczoP
zxmK72QNz5LHtPlfPBeb/8V7YfO/sEs2/zuI5lmb/8V5d5rnDPG8Ppb/fYn8
Upv/7UbrApv/vZPWZTb/+wqti23+977o/vmx/O+wCLf7KDN5eM4yn5pnyWPa
/O8edLX5X483xPK/e6jvUWHwulj+1+Npw3OWeN4Qy/96vy4Vy//Cn7T5X48n
Y/lf79cti+V/vV+XiOV/vV831fCczcMzeNwWy/96v2hLLP/rx0ljLP/r262P
5X/9/2+O5X8LyZ5oPpuIh9pY/hfzps3/ejuzJpb/xfrI5n8Hk/9j87+X0jxo
879X0TpR89lk+JT87wayJzb/i7ilzf/OpvnU5n8Rt7T531HkF9n8L+KWNv+L
+KTN/yI+afO/iE/a/C/ikzb/6/20wbH872kUf7b53z9T3tDmf3Fep83//pby
jzb/i32CNv97GJ0XbPO/OH/c5n970PnvNv+7NcrDjorlfy+I+jMslv/9dYRf
Ecv/4nxSm//9ivbr2fwvzue1+d+z6Xxk/9gk//sWnWtv87/VUZ5rXCz/OzDq
z4ux/G9FdP+gWP4X59va/C/OF/btSP7X0fnOXci/AP485fFt/nck6QFs/ndN
1J9HY/lfnHds8784b9rmfx/zf4jlf/3+wWdi+V/kqX3/Jf+bpfFg87/e7KRj
+V//eRtj+d8utC7oouwA1inr8+z/XRvL//rxtzqW/y0m+1io7ADwT2L5X3/f
4jz530Qs/+vf2w9j+d+PyI7Z/O/1ZB9s/hfnsVYr+5Dk83CLlH1I8nnEXZR9
SPJ50CXKPsh+W5v/xXnohco+JN3+lK/U81LS+fFQFsv/3kh2zOZ/f0T7N23+
9yiyG/53yf/iXOASZTdSfC5zkbIbKT4X2+Z/cS65zf8uoXPhtd1Imf2qkv+9
OXrf58fyv+dG7/XUWP63H9kTm//FedMl2p7wed82/4vz1m3+9zY6797YEwc7
Y/O/0BXY/O9jZGds/hf82/wvzoW3+d9oG3HvDbH8r7fns2L534fJ/tj8L8a/
zf/6+hW1sfxvmvYL2/wv7I/N/+5S58ba/C/nbbkd5G1t/vcOdV6tzf/yPmL2
iyooj3ywr0fNn9th8us1od2h+tU8X7/3xRUnh3ltqnfN8+maxYv/GebBqT42
46fuf/37Yd6c6mnz/Lij9vpuYZ6d6m8z/temf98d5uWpXje/fyc/sP+IMI9P
9b0Z7/jyGd3u3n816oHL+zr19gWhToDqh/N81/LK591DXQHVG2d8zvNXLwx1
CFSfnPHBw0uvDXULVM+c57tuy24fFeocqP4542PmnJYIdRFUL53xTld9ty7U
UVB9dcb36PVijxWtdozqsfP8WHvkYeeHOg2q3y54y++ODXUdVO+d8ZrN97cL
dSBUH57xDUsnbwp1I1RPnufTDz9oKA11JlR/nvGZ9138YKhLoXr1jN/66PV3
hzoWqm/P+M0HdWkOdS9UD5/xm8ZtGxDqZKh+Ps/XzxQtvSPU1VC9fcaTZy4Z
FepwqD4/44ceuHheqNuhev6Mzyx5ammo86H6/4y7K/9+eKgLovMC2B/YePIj
+5RfU4nzBRi/fdWDV4e6IzqPgPGt1Z9fGuqU6PwCxluO2LEj9K/ovAPGj5tS
8bPbpgQ4H4HxP3X7+IbpvXGeoJy32C46JxTnLW5l/LkLKzqEdp/Oa2D81HZn
HhDquOh8B8YLo/MBq3EeBON/brr901AnRudHMI5zLem8Ccz7wYRLTrk71KHR
+RSMD378xXdD3RqdZ8H4uDsSA8K4GZ1/wfiKK3eeFOZB6LwMxj979MyGaj73
tpHxS654bcCsohk4j4Px313WXB+eQ0Hnd7A/8O2v+s4MdX103gfP1513/XZC
+FzofBCeTy/68ZKOoW6QzhNhfGTNonPCcU7nj/D82HfwgzeFukQ6r4TxFUd/
3j/UMdL5JjzfXX/wmLPD95rOQ2H8+76fzAh1knR+CuPjJv23W+iv0nkrPN+9
fF5mZaTD9OezML5+/YUzQ7+FznNh/KjB774WrrMO8ue/8Hx3/ek1z4e6UDov
hvHDTj63W6gjpfNlGF913KTSUHdK59EwPvzoAy8PeaDza3h+rN337GtDXSud
d8P4n77Z78lQB0vn4zC+eFvl9FA3S+fpMD547m+DUGdL5+/wfHrUsdPODfWE
dF4P44sG1n0Y6njpfB/G+8895PhQ90vnATG+z94XvB7qhOn8IMbfu+zeg0Jd
MZ03xPP1Ubd1qAh1yBf484kYv7ngvtmhbpnOM2L83dHr3gjXRz/x5x8xvvPE
M/8e6qLpvCRp59mz1ob+DJ27xPPy8KU3/C6Ms9E5TYwff8jzc8O4BJ3rxPjK
folfhXEJOgeK8YHvfTY+9Ofp3CjGf7rrmA5hXILOmeJ5f3KnA74M9ed0LhXj
O6aF/sFYnGPFeI/C5IIwLkHnXjE+JH3qO2Fcgs7JYvyBnSf8LDyfi87VYrzd
z5u3hHGJZtr/QTjpbVA/UvyNBl13hf0K+F3YFwIcfpr3M8Xf+Dv5deQHMg4/
kPxG9h/gN+r9+pkAfub3an9JJtB+qfgbY8mPJb+Xcfi95Ccz3pv8ZOyvBX4r
+dXYtwocfjj57YzDbyc/n3H4+X4dJP7GA7QuoHUE41hH0LqD8f607qB1CuNY
p+j9bSmqD7YW6yDGsQ6idRPj19G6idZZjGOdRbwzDv79PgzxN/rSOo7WfYzv
T+s+WicyjngXrSsZx7qS1qHSDq1Dsd8aONattM5lHOtcWhczjnVxkcoPJAOs
o2ndzTjW3dhXBBzjf71aLySCVbSupzgA46UUB6C4AeODKG5AcQbG6Yq4hPgh
FJegOAbjiGPoeiHAK0wd5QTXUaS4CuPVKq4i/oa/Ip4g/gbuozgP4xdRnMfz
ts18LurdgreCAv+5I85cr+LwBQWj6XmXqDg8dAXTiR/RP0DPWa3i8AUFxRQ3
L1ZxeOzvXkJ4mWlnOXCn26ngeBJw6BmqVbxd9o8Xq3h7QQHiLCUq3g4eqolP
yftjvB6r7HbSfUbxN3xv4K9TvO53yj6n3GsU39P2OeUQD4QuDzjih1qHmHKI
Nxo77L6m+KSxww7xTGOHHeKfxdoOcxwDzwE44qtaF5B2iMe+oOxtxu0/3Mdv
Byh7m3FFZIcvUvY249ZTfFjb2wzpzCci/sx20VH8+RJlbzMO8WrUCwSO+Db0
O8ARD5+p7GrWPUrx8xeUXc26s8neDlB2NesQj71I2dWsO5ni+Rcqu5p1iP9D
VwIc+YJ5yq5m3STKLyAvDfwMykfcqexq1i2n/MUGZVezXLdwprKTUnf3RWUn
Ued2KPIyjK9TeZn15v5i7L9kHPsHpys7KfV1kSfX7Uyg+K34af7f32I9i25/
IusFgPv3/V3si2Qc9cG+VXayKbiF4kvWf3ue+mP9txS9d9Z/eyAaWNNi/ttB
5BdZ/21bNH6eiPlvvi7f6Jj/5u3M3Jj/5vd3zI/5b8fTvh7rv51P49n6bxi3
1n9D3M/6b5f7jW0x/62B6l5a/20I1d21/hv8Deu/fUX8WP/tiOh9GR/z3/y+
mMkx/w37Ga3/Fr0WvRbH/De/b3FZzH/z7+8nMf/N19vk58F4B9rXY/239sP9
+LH+2wDaf2r9txTtJ7X+2z9pP6n137A/1Ppv3n4ui/lvnahepfXffPx2Vcx/
8/Hw8pj/VkX72qz/hn1G1n9bmcd/g87c+m/YL2D9N+wfsf6b9xMSMf8N+yOs
/wbduPXffLvJmP8GvaX13zy+Oua/+XaqYv6bx1H/Y5tpZ33Mf/N2qTjmv43J
47+9/H/If4NO/v+W/wY/4/+W/3YXvV/Wf3uW6kVb/60z2Rnrv50a2cmZMf/N
r19SMf9tZWSIq2P+24Vkb63/1pvmd+u/XUn21vpv3v6si/lvWGdZ/83vN0zH
/LcCb7hj/ltP+r7Wf7sn+twFMf/Nf+7HMf/N+yHVMf/N6y42xvw3vz7dFPPf
/P7QzTH/7bbIvr0Q89+OoPnX+m8jIj+5LOa/PUN1fa3/5vdFlsf8N8QrrP/m
/dh0zH/zdrsu5r95PBPz33w7DTH/rSG6b0bMf3uF3kPrv3kdzicx/w11ua3/
Bj2q9d/896mN+W++35ti/pu/f3PMf/PPuz7mv/nfG2P+m79ujflvD9L+UOu/
efswNua/DaX6ANZ/25viJ9Z/20z+j/XfplFdZeu/IV5k/TfUbbb+G+o8W/9N
14UW/w318K3/1lHN4+K/PUnj2fpv15GfYP031KOw/hv8WOu/oc629d+uobrc
1n9DHW/rv71Jdb+t//Ys1Qm3/hvqilv/DXXIrf82nOqWW/+tE9U5J94Zv474
t/5bmuJd1n9D3XXrv6FOu/XfUNfd+m+oA2/9twVUN976b/C3rf+2P9Wlt/4b
1oPWf4OezfpvGP/Wf/PtVsf8N4+vi/lvqNtv/Tdd53+9wStj/hvslfXfMP9Z
/83jS2P+mz7XQPw3nINg/bdZ9D5b/+2uCH4/5r9tIvtj/Te8P9Z/w3Oy/hv8
HOu/oS6c9d90nQ3x3+B3Wf9N+2niv1Urv078t1k0/qz/pv3GOtPOvJj/9v/0
dq2xehVVtBIRKAYI4g+CUEhMFFRiLFHA0FFLFH6IJkIElaJifCQqWN5IRfBR
SCgUgrSEYngZShVogT6gwJyWcml7C57etpdDSy+n7bnlUD7lISWmRPH77pm9
F2ud7/5T509zd07n23ufmT1r1t5zxs6fKn4zHlLx2wnp/gvFbxem+zIUv9n9
Gorf7D4OxW92f4fiN7vvQ/HbQLofRPGb3Sei+M3uH1H8NivhT8Vvdr+J4je7
R0DxW57uT1H8dmi6b0Xxm93PovjN7nNR/Gb3vyh+s/tiFL8NpPtlFL/ZfTSK
385O99cofrP7bhS/TUv34yh+Ozt9f0Pxm92/o/jN7utR/Gb3+yh+W56+E6L4
ze4PUvx2ceI9FL81fn+1hd/sniPFb824ebmF39I62sJvzb87WvjN7mlS/Gbx
RPFbM283tvCb3Rul+K35udUt/PZuwmeK375O9zv4Pehx+7u3HPbqVKtrrFx+
0zHf6Ew91+reSpd/6YwDn5h3Nc5Rm/zNKwe7A9rW9dzldy2Y2R0SiO+uz6ap
3feF+/dMfuhpWz83PNPmYeHyy+4+44XePE36uHzhG7ftnj7Rvrtbufy8edcd
37XYx4fJv339nrd7uJTtreNxB7/efXyn2FvH2QdkXTwwKvbWsTNxdjcO7BB7
63jqPt/tDultYm8dTxmTvyR21XFpNuXBA2etF7vquHX6E6Ojh6wTu+r41R/8
66c93Mt21fHyP19edyGEjwOTP/KR3X/r4W22q4ojz5zfHRE7xK4qvvyH07vx
f5vYVcXXzj3+PfnTzOWvj8m3iv5V3O+tkU09PMz6V/GYT/78zh6OYv2reOvq
PS/0cCPrX8VHf/GZ7WdOWiB6lnH7kYeO7ftYzzJ+6qqjuoBvRPQs45zsoO57
3CL6lHHfMT1Xij5lXDZ6+eDg0ctEnzLGgfOGB4+eL79bxCtOHOpuGfT73UWc
cvgrm2YOFNJ/EZ/f653Hpxy3RPov4qrD7l/byx9xP3kcvWzFZMwL1KuvPn/X
k4um2DlZW8ezcOM1u2aj/hvfi3n1EwdP7q0v5l+TT/3rCZt69VFmp8nnTf/e
JSfsvUn8n4XdH762u+RtED9k4WuPLuyG7vWifxZOG5Pbfs/mexZW/+PqH2E8
Fy7/5Zg/7b4AH3fhrL2mlhjnlcsfv2Fbd6FdLH7OwpQTJ+zA+O+4fL8Fvfky
T/yTh+FvfvSHPR6V/ZOHa/54wYIDZ20Q/3Tx0ujKv0+fOCT+yZOfc/FPHk4c
e/458U8ebh/z86D4IQ+dV055cMpxT4of8vDQ+3+fYXx63A2XTnpqfW+fy37I
w4wLJwxg3Pr4CR86fXc+ePRcsbcIgz+7/8oer872FuEvD56y7sxJz4m9RVjx
m8/eevWcQbG3CMWYP9eIvUXojPlnQOwqwv77rLmjt/9lu4rw42mnXjRx2gM6
j8LjD0y47z3j3+WH//u8RxdNuVn0L8Pozfsv7+VBWP8y2TUg+pfhgrW/29qr
c2P9y/QeV4qeZdg244NfOXfGfNGzDHvqCdf1eE7Wswwfn7/rzjnfmK1xOMzd
8q13j703E32qcOTlz92H8eBxOBx00sun9+rZ+HercP6X9xrZk8+V363C9+9Z
e8g7+TXSfx3OOufGufC/r1/hgTH/PCT91OGif8bRgT2/luc74caTT58J/Gd2
dHEFxSXwRYxD7Pk6crzKXc7z1PqpBLeAL5pM67v1X0WOe/5dgMjzHfGS54X9
bhm/SLiodjnjB/BFvP6aPmXkeIv5J3HG5ffTfCxdvpPGuelZxDcIp/l+TfCM
xwHBA+CLttH6a/oX8W1eF1zO8bBw+UqKGx4HIs9H8EVzaF74vlVwpsfzhLu2
iV15wi2l2JULfvB4ruu7zxte7xDvOc77+I3PUzwEfuE4Axw9iea7x/PI89TX
fcHVvu4L/vR1X/Gby28hvOTrfjyJ8Iyv+4JP/FxZHGfdj/3Xuyz2Xxey2D/e
ZrF/3MvSvNB4lcWbKC6B9+B9h7/fsJbwic+jwPsR8B5LCM/7vBA8g3nA67vH
k8D7Go8n4UXeF7h8X8LbHk8EL3k8Efzg8UTWX4+rQfZTLuf9SOVyxvngPRhv
exwW/IbvszCe8TgseAB1nSWtv76+hO/Qvq92Oe+bfD8YeD/i+8HA+wLfDwbG
875+hZMIf/r6JXgM+4PlhGfwXQ/GDx6/ZX33fX3gfZ/v6wPvp3xfH3hf4/v6
sJL2I76vD2to3+HxOEwkXO3rpOBM39eHIxinufxjhIt8Xy+4Bd+JuJhxiL//
hhcDX2B+ejPVsfDzdm79JdEf9eHcT5XyEfeJ/03+qI7zVA+zVfqvUr1uKf5B
fTj/bpnyL/dqfEjPL9X4kOqCVog+do57s+hTpnr1raKPfQ/c102XW30461nE
y8Z+908aJ5M+i0XPIjb34UI/fv5p0b+IzfeNFR8W0fLgrH+R7q0eEf0Lsatd
H8525bH5jvTdYpfJHxa78uSHJ8Que/4p/z2TNzz8arE3j833rhWv5um7AYXY
a3XpW8TeXN4v8AaPQ+CN/uM/S/3cKX7I0nhYJH7IUn2XxlV7foX4IUv1eLov
y7xOjP2TpfzpRvFPlr4noLxiJuMHeIPnBfAGz1/gjWkUN4A3lvSNP9buEL9Z
Wyh+s7ZM/GYtE79ZWyV+s7zYM64fP79O/GltSPxpbaP409qw+NNaIf60tln8
ae1F8ae1EfGnfQdgjs6XwOME+G1nqitOvyvyu6WfIsg4dDnPU+C3pk7vBo3P
SW73GAK/bUv5X/5dO1+vuLEMEjdcznES+K3pf6boU6Xz5nNEnyqdx7d5Xcrz
qmeVvuv+mOhp30V/QvS0enLYyc+vEP3r5LcrRf9a/Az8xnaV8rz5H/it6vve
65SvXCR21YHjPPAbr2uIo7Iuu5xxAs7jf75/viZsID8A171E7xe4biP5B7iu
ovEPXDfSdx51/PsG7IdO4HUQuO5SWveB694i3AJcJzjK5a8RTkPeiscV8la3
0zxF3qp5j8+qH1L8H1I/yDqLvJWdO5tAzd7biPon9b9d/ZPqQCr1T/rdneqf
OFZONa1W/yR9NH9Xx+b8wirxQx3H0sWvDoofDN+uFz/UaR3cKH6ok54a9+39
vCj21glX7BB7rf9K7K2T33aKvfa7GB8mZz8gz9WMkzVibxWbegDl56vkhw1i
b5X0GRZ7zd+bxS7rZ4fYVaXzjKNiVyX2Iv81l9478l+XpDw462+4XfMvpbxH
a+D1Wc9S9EdejN8X8mI8bpEXa+ZjLvoUaTxsEH3MvmH53ULGD/JlPI+QL+Px
bM302ij95DJPreF+3vR8+jsLdk9E0tPlvK9Efs3OE5m/TN7Ez+Vuj8kb/VeI
P7NUb6P5iCzdg7Na9U9/G65D3u0c8g/ybvxewO/y+EHejecF8m4cH5B34/iJ
vBvjH+TdmnVhkfgtl/UUeTdeF5B3+xXhKOTd2J/WzE+rxT/2nahc/GP9rxf/
2PevNoh/cvEn8nG8T0E+zupU2Q9FaMrVHxY/FGFyqi9lPxThEvIP8nE83qyZ
3avEXrsnqHUeJ3BcRZ6O5yPydMfS+os8nd3vxnaVYi/ydDwekKfj+WVN8xZY
hy+h9eI9+TtaN5G/u6IvrqjCmfS+kL/j8Yz83QkUH6xVqZ5nqehjz68Sfaye
X/fLte4vXM74ypr9zkLpvw6X9R0/Vu90m/SDc4jcT6c533RUR94DzgPa3+l3
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

   "], {{
     {EdgeForm[None], GraphicsGroup3DBox[{Polygon3DBox[CompressedData["
1:eJxFmWXYFVUXhkffOcaMetxiKwoGqSAqKggoiAGKEqICJi1KSaeoKHYHKCqN
gIQ0SHdLh43SEhZ2fPfjM9f1/Vjvvs9+95kzs2OtZ60p2rht3TZHRlHUkz8F
tMOSKAq0U+IoqgCsoX06jaKr4NVwP7ga/AX8CnwlvAp+Ci4Lz4N7wRfDc+Ae
cBl4LtwTrgp/Dr8MXwd/A78KXw/vgF+Hr4Y3wM/C18Cfwi/CFeFP4GfgKvBm
+AX4Wvgz+CX4Cngl/CRcHf4Wfg0uBc+Eu8IXwbPh7vAl8Hy4N1wOXgT3gRvA
RXNR9AF8D1wM/hC+Fy4Oj4Xvg0vA4+D74ZLwePhO+Fx4BNwYvhj+CG4Cl4En
wpXhTfzW83BD+Dz6R8GN4PPh0fCl8GLGPAaXh5fBfeHL4aXwE/Bl8BL4cfgG
eC/8Bnw3fCHXGQPfBH9H/1vwjfA++E34AbgUYybAleCN9D8Hn4o9Efn+a/P/
6fCL/K8+PBt+Da6pZ4KfhuvAM+CX4N/1P7gk7W9YDbgE7QKsCBzR3sb4aeAL
8K3wVPh5+A54Dvw6fA+8BH4HvhdeCg+E74OXwe/C98PL4ffgBvB8+C04guvS
XALfCc+F34DvgufBb8L/6r7hsrHn5mr4rNhzUxE+M/bcF4K/K/C+OAneCy+M
fUZ2F3juT4YPFHhNToEPwSfwWw/AVWLvxzPgXwu8hqfCP8AnM7+PwYVzXrcK
8Bkay+fWcA34B6waXJT2e6wqXETXYEwr+Ab4NPgh+Eb4dPhh+Cb4aLgBfCV8
JtwGrgmfBbeFb4bvhhfDb2t+eKgifB7JfdUI/u3+qc9racYcRzsD/3Ca5jDv
s3UZXCj2+SsFp7HPa0k4oV2PFYWPpJ3Jd0/XGuV91i+HT6b9hN/4JvM5xfjd
qbQdU/ubs+j/s8DnvjD8LzyG6xyjfZK3H7oIPp52O1YODrH9x3lwAW1xrjmN
tlNqf3Mp/SfFvpau2Y1nnsQ1T9Q+4ZoTMx+odZ4dvEfX0T8neO+uh2fBNeC1
8Ez4RngN/DF8E/wJ3AK+XPub3/0IvkBzrevDpbW/4bNp2/Fbt2hu6a/G59X0
T4Iv1b6Hm8JldY9cpxl8CTwZbg6Xg6fAU+AKOg+MnwpXhFfA0+Br4ZXwdLgq
vAqeDF+mswTfAZ8DD+c6H2vfcD/fMzcbEs+n5rV08P+6MaYkPAPuAtcN3q+D
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
qWs7n9LfMHUd6XO4TurYtwWumjrXUL7cIHVd6LO8tZe06N1cp27q3GRr3ppc
OvNLuFHqutMXeetP6cNNcKnUWl257TU567G2XKci/W9FzoVr5az9etFfIXX8
VV58c866tAf91XPWex3lU3PWTp0U01PrQ8V51fNUK7tXcTy1/lfefUvOGrgn
/ZVy1loPw6UZ83TkvP7q1NpA+f75qZ/3okz7SXMql6+ZaYnuigM568DO2m+p
62C78q5Hal/dqdiXery+d4j9cAs8ljHfJ46t4/OuH0gT3s/4vxPnI/LhvyeO
3VPhKjlr2jaM+SOxJpxG/w+J4/4E1X+yvS1fujtx3BlF/7+Jcz3Fi32JY+WH
eccC+eEP4L2JY9aYvOdR8/mQYmvieDea/oOJtcG4vOuvOnd3MOYf+hvB0+n/
NbGumAzXTq3BNuddu1W+oxxKNS3lmMp9aqbWYBsYUwMeEbm2o1ppd7i+Yntq
jac6z2+JNcwU7Y3Uem9j3jVa5XTK3X5JrH8m5e1LFe8m5l3TVR6nvO/HxJrn
o7xrOX3h+2LXhpUn3iMdlcUm+eRWqWubuxl/XOp7U73irWCNpXl7Pjjeq+aw
LLg2qvO1NLhmqnP0QrCGUI3i5WC9pfrMc8GaRrWIF4M1hGoULwXrG9UfFgbX
T3VOlwfXW3XWng3WT6qBPBmssVSLeDpYk6m2szpYpyq+PxOs21R/6BecS6gu
tCJY4+osPxWcb6j28kqwFhwOrwnWTNIMXyXWfuovlPkf1XxmZHUb1canwxdG
ro0r1so3qoYgzaJYfx7toqyu9V+dBL4gcv18SXDdWf5qTeKasN5TTElcS9R7
ganw+ZFr8qrf67s52tmJ6596D7IouNYsfzUnce1U701eDdbNijVzE9ct9W5l
Hnxx5PcaqjmqhnME7cPBvkM58qzENVW913g9WHPrHK0Kzg2ked4NrmVoj/UN
1rWqCw0KrhfojHQM9nHKqTsH+0Hl3Y8GayDVNJ4I1sGq0z4WXLtXXbRPcIxX
3bVVsL5RTvp4sLZWTbhdsP5WDt4m2I8ov+4Q7L+Uv7cOjk+qD/QOjp2qz7QN
9tfK2R8J9rnK6x8M1kzKWx8KjuXKi9sH+2jVARYH17IVs/oH57HyOb3gxjnX
cAYH16d0locE11bk64YG17Dk3wYG16fk31YG50vSVO8E17PkP98LrrPojL8f
XH/R2e8WHKdVA+keHKeV43cN1gqqmbwRnAvJT/YMjv2q57wWnBep3vJmcO4k
HzgguGYkH94jWDeoDvN2cE1NvrFTcLxRzaRLsI5UTWZ14ncZei/WMlizqnah
92Y6I0fRLk1ch9d7qGVZDVPv1OYFv9eRr5sf/O5HsXVu8Dsh+c8liev/eme0
MTg/30P/+uAcUrFpVeL3Gnp/Nya4Bie9uiE4t5R/m5zVLbX3PwyuCUr3at1U
a3pVez24ricNMyq4fid92zp1jrwXHhdcO5NWHxtc85JOHhZcK1S8aJc6d94P
jw+uzUlXfxKcg0krLgh+F6U4olrfysjPJb0o3ah3heuCc3Vp1NHB9Upp7G3B
eZf0p96nyVfEsd+jSt/+jTZenPgdhN5vbg3OzXQ/m4PzyX26H8ZUivyOtU3q
nF39W4LzN2ndtcE5ubSx3jeqtqk66qbg/FNz8j+uiBzc
         "]], 
        Polygon3DBox[CompressedData["
1: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
6rlH6+SvsVXqk8fq5MuXq+ghazxrPWu6rZrUoNaiauL769QY1hq+o+9qTWVt
5TPurJNzzb1+xx11NLHa2Brsrjo5w9zRq40f5lBzqb9xdB2/wtrPGtD635rH
2keNfXcdv0YtMaYTv8X4qL78sNU3xmC1wP9jKG1ncM+CnP/aI/PLHKCXYgwy
Nxjz9WJ6t/nGmtba1pp3a/uX423ODy2JH+Yna+NF2vn/Sh39+G4V/WR8VPuo
gV5vUtNb21vTP1BnDKZtv8+xsQa0FnyoE39Gj0CvwBp/kyZz1lz6Ridz2Zra
XG3O3sPas8rv932MrdbY1trO8T2b1KTWpq6BHZvUNNY21oz31anZrd31GPSD
nANqY2tY54Yeg16DNWv/JjWwtfDznfhlanq1vTnH3KOnobehp/FgnZrB2sEx
errOnOvR9p9zUc2t9l4MfryO5lZ726dP1OlTtaka2b6eUKce+aiKPjXH9my/
39zrnDI3PVZlrrnG1TofdLL2nWM9/K4qc885WvldVeauc9gw9GSVuW0f1eay
Kn3nGmjMZVXWhmtgBviZKmvDOT4d/FSVue+YzwiPqTIXHK+p7Xx3/I0Retcz
tOvJPnCtvFKlb8yRU9r14vhac1l72QfP1Olzx/K1KmNhH7tWJ1Tpe/vQtfNq
lb51TGZy7lcZK2O+ff9OlVzgGpvesayy9szhjuWkKrldzTRLG//NP+ajXm2+
MN9MqpPrP6lST+q56b2pkd6vkyPUBh9XyR3GUL1iNZix1Zzj3Hu7Si5yjTlX
3qqSe3xn19LEKn1hTrO++bZKrjPnmou+q5KL1Rz6ba9X0SL+Rn/rl1VypTnW
3PdVldyr5vq5jW9qsXPaeDdPO78vaOe789nYo+bUa/y/xqyTw52rX1fJ7e5B
6KeUEj/k4BItrCber0kfWXur0ey7AW39sV3J2uzd1o/LlGjPrdp6ZfuSWK2n
qv/1mxKvVc2j9lmsxJ81J5mbli3xr/Vw9Qu3KamP9Witn7cu8W79fmPBtiX1
oTHWWLtLiR+n5239tVOJf9e/rQ8HlMQKY5hadsuS2GYMG2kuKYltetzWdzuW
+GsLt/X00iVaW8/b+nDnEn9u67Z+26EkNqvx1P5LlWg/NdgIeO2SvQE9d73k
tUq8eGuSUW1/H9WkRrJWsiY6oonGtRZZsET7GlOsfZ1Txho18rGcL1yindXA
x8ALlWhjNevx8CIlWlYNfBzcu0Qbu4ehH7FBiT/vnora8B9V6uUpdbTfn6v4
H+7h6FcYhNzvUiOoFVYt8UvcM9HvX6XEX7GmGt3mB/cD1F/W49OX+Pc/1NFe
/6ziN7h/YD6xRnuryZ6dWvOnKv6Imkfts2FJPWaOUQvsU5J73IPR21+/xI9X
cw2BVyzRYnp0enXuITRNNKveoJ6SWlbNptZfqUTLqYH1zvUg1cZqGmuN5Uu0
jhpQr1OP3cCrhtarHNiJtlbzqOVXLtFCajhrmxVKtJ2aXS24SYmW79PW131L
ai9rMLXfuiW1mTWV2rJfSa1lzaa2Wq+kllPjq8X6l2h/NZm14RqtVrP/1Hob
lWhJNaPabuMSLTmKPviqE42nPtBj1WtVM9dNNLzemx6v2t6crRbYvSSXm3PN
xbuV5GJzvFpl35Lc7x6K/vqmJfsh7mHod29Wsrcxrs33jqljq+ei9+Ka3LuJ
5lJ7bV6yv6ImtbZavUSrOh+tjVYr0bJ6tHq1eso9m3jEesX7dFKPWp+6d6vn
WzXxfPV+rQFKE41uLb1miXbXA9YLtkaZronHrNdsDTJ9E89a79qYOkMTTT2U
8z4lWtua7gR40bbW06PSq9KjOLJJjWetPH9J7ecetnPrP1X8dveg7ft/V/H3
3QN37v6ryn6Je8bO1V+q+LXuf+sd6yG7n+v61u/+exV/2TWl97htJ2vNPV79
3j9U8W/NifrLU6rkSteo3qBr2LVrjNBb85k++3B+/0fGshJ/R02jthlc4t8c
DL8LH1biR1lDWzsvWVJbHwq/b+wr8aPcU9av/1sVv9w9AvcK9Bw+rLOHrd/+
1yp+tTHzY2NNSSw1ZundGrOMZX5G/31qlWcZw/q188HYZkw4kfM5SmLFdG08
naVk7Esbj7slc6VHO14zlcwt97v7tuvf/XE1/9zt7zUfu6ftXKxK9gesGfXb
f1+llqzb+TF7ydrr2cbrWUvmsnrCWvmbKvsnW7fxzPVpbHXP3/2ZH6v4+465
Y//HKv68fyPg/tAPVfYHpm/zR6+SuTxDmy9mK5nLvoPv8t8q+w/G8M/h40ti
+yHwe/DhJX6iOewz+LiS3Kanp7enJzu4yZzQ+9MTdK7oYR3K+Zwl3paeq97r
GGNGE0/22lbPq0f02AZxPleJ93aQORMeVOK/6rnqvZqj9zEec/6D/dtJjaTn
ph8xd4kX91+00/fO/U5qfOsxz6dwTG6iV/T2Fi/xbvRc9O6WKPFi9ID1gtUc
A5v4OX7XVI43m3iA+inzlniDekL6LfOVeEV6RvoJ85R4ScfAP3ZSI5v/1GPD
2/rI+sq/EXE/6S9V9nt8Z71nPV77whjyKefHlsQWY6i12MCS2GqfGmsOLOlr
/+ZA/3fvkv1Z/8bG/aCfq+zXuOesX7xXyf6l76g22b/k3R3TT+BjSsbaPnes
DygZC9e43rPPcO27R+1v3bNkP8QxVfscVDLWxohb2vrQ2OGc0qv2NzrX/JuI
seaGkv1IY4Tete9o7DBHOLf2K8kd/wMitqDW
         "]], 
        Polygon3DBox[CompressedData["
1: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
l11gfwM8Wtz4
         "]]}]}, 
     {EdgeForm[None], GraphicsGroup3DBox[{Polygon3DBox[CompressedData["
1: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
t5e/jPYjB+g/Ptqf7oefifa5BcbzdLTPLQefAd9Jn2/pc3q0nvgGPis69n8P
b0isS1ZqPWnfy89kl3gJ12rtV65l2s8Z/yVfkPX5IWv199LsflsSa51NWX/d
Y3Pm9/bk7PvE0kkr8bEPaZ3gY1dVdF63G14O3yffWex76NnLE2u9n2gHptah
v8CnpNaMP8Mnp36+xlEmezexdFf70tYt0mZrsjGtyMa/MXuXNZlNY3uVOfpO
a6a0deWf8NmptfBv8GC4YpHnqQJt5Pow49r0+REeoBjF31Xhj4LfZXk2nyt4
xwf5/53F1tq/YhuUWoPPhU9PPfYN2btobOuzOVyXza++hXTufNqztOe55yyt
nQoet8bfn7/bZeNvW9r69w/4jNT6dx58JryM3w7XPi32nCzN1sBi7Ldpf2G/
KlrzKS8eEq355GOvz3K0xfANWY62RL4ZngrPg6+J1oIL4aHR+m8B3Ddat70P
Xxedvy+Sb4/WZO/CN2b54FL4juiccTV8GTwB/kP+PloTfA6fGa19v4Nvj843
V8FXRGvN+fAt0fnpCvjm6NxzOXxTdF65DL4tOm9dCf+P9kTm46e8c69P4B6p
c7KP4G6KB8zV49qbzNVGeDKc8C02wE/COXg9PAk+SJ8t8Ez5Iexr4XHyD9g3
wzPkq5TPwrfLt2NfAz8mX13snO9jntudv9dhnyj/gH0BfKV8eLFzzc/oc3xq
Lfoe3CV1XvgpfJzmmv43aO3Qv2a2Zk5KnXd+A/dNnTt+C/dLnct+BZ+YOvf9
Gu6T+jsfDX+ct++pAT8FDwv2k/KXeqeT4N+wTwr2b8qFx2Z9NMZHg32v8mK9
d2P4RfpPCPbbynknBvtt5bZPwC8lzovHB/t55cJas43gF/jtI8E+XHnxA8Hx
S7m83rsO/LTyBexjEseaK4M1m2LWg8H+Wfl+/WANrHU7huvZxDn+jGhNL505
Upogcd1A2vMtPYvnXBWsaZVnjQiOp6oV3KJ4mTiPrlRk36C9O1g5vHKDLB+a
gP0w2mej9bryROnuLtg/U67Bb8vCcwu+z/fwa4oRwetEOfXU4PGoFqFc80u4
d2q9rHdpSnt38LNUD7k2OK4pJl8jbZE4r785OI4rZ38yWtNLb98aHJdVV7k3
OI7ruyuufg73hO8JjvWqsTTE/gXcC/slwTmCNMzQ4LxA+mFIsJbeSdtM2iix
D9nBezZNHONKcR0Fj5LvC16H0lH/Cc53pLtUSztH/jfvPdcf/l25M/83NXE8
vZ1rZuI4e2NwvqB4nQbnp9J+B3lu+8TarF5wfiTfcLp0YWJN2C7TpdJrpYNz
XunAffy2TWLtuhNunjgW74JbJNa3B+B2ifXhfrhtYv15XTaf+hbXB+cv0l03
BOcm0m9nB9crpPeuDs6zpJf+G5ynSMvJN1SHJ/Pu5wTncdKEFwbnQdJO5wbn
ldKHFwfnRNJXFwXnPtKW5wXnp9KKF2TrRDpI61a1i2WM/a7gPFf1wDuD81bl
L3cE57nKR+4Lzm2l8++Xn02s/4cH58XKL84KzrWVy6hGoXrIFMY/JfV3/0j1
h2CfrNrUIcF1A/l51U+a6BvRZ0awD5SfqwUfn7j2WxM+IXH9rVK2lnTvsnCr
xD65NtwzcZ22DtwrcZ1W9dd+Wm/c/1Ds3RLX7uSnj4U/x14de9fEcbkc3Dpx
jLiYvfNo6ty5AvYjE/vwYrhD4hgxJYs7P3KfMthbJvb/FeGOiWOK4onG9p1q
I8G1ccUdxQS917fKi7H3TlxzVh33LHgp9iuj6zKqIVcLHrPi0ZnwtYnr5NIp
Gs979F+U2g8U0V4WnLupFj0oONdTvX1AcM6ouu7Jwfmg6t6jsu/1Ifc5NThn
VE2+arAfU7y7HL4ncR1bsVTz+SX9+2K/InENuW7wnGjfqQ6kWvSBvLWk5uR9
5dHBuarqzNI1+o5vYD8/uL/WquKtvssX2EdE176b0//q6NqTat2dtH8Sny8c
FZyPq878fLT/0nnBi8G+VLXT/wXHR9UYXwqOj6oxjo2u6aiuJV2gb/FD3rpB
+qE192pJ/1MT10VbyX8lrqm2hk9LXDt9OTgun5JpPMUL6bppwbmDasVPBecg
qksrLkjPqt6g84058Gnw9ODcRDXh2cE6QTXe54K1gWq/rwTHfdVUpce1bt9k
zC2wD0p8xtEGPj1xjbc5PDjxeUfX4FqEavW9gnNqnQUcE1y7UD3/hOA6gGr+
XeRHEp/ddAj21arbtw9eq6rbd5RvSnwOcqT8XeJzh+7BtQudFxwXXGfQucCt
0bUBrZGGwdpD5ymq72otbc5bl2kvfw83CN7LOtNpFBwjdO7TBB6Y+Hzn2uha
mOrDbTP/Lz/fGB6Q+Gzo6eBvqrp0n+DzFOnMMdG1M8WyxzJfJM2zvWAfJV+1
Ca6XuEahupByF+mA5Vk9Wf5vdcF1XdVw1sEVE+s3aX+NYQ7vsqJgPSAdqDxF
66uDagMF15CVd48OrvupzqM6R134OX47KLruplri4oJry6qN946uh6qGqRqG
xjmb/isL9ska25KC69LS8//CzRLrcNVCFO9ezbtOWRuelfc5mNbPQnhWsG5R
zr6Z39ZPXPfIBWs55YkJ3DBx3jQ5WFco51UNRrG7IffbWPC7qOaztuC6t+oS
qge8D3dNXftRbH2H524ouF6t+o/O9LSuluXtmxWXX8v7DFBrbwW8hf4NEufU
M4M1s2KK6rv6FmPyzksUx9+FX6DP74lrCzrb0Vqam3e9SnrgbXhNwTV51QR0
Hqi9/xf2ksFrQ77z0Gj/Kx91fPBZnnKiHsFnfMq5FhVcf1beVDnzk4ohPYPP
+5TX1MxygWNT15yko17mWacE136VWw0Org+rNjUqi6Ev0WdgcD1Z+VT/4LM2
5Uejgmu/qn1VCa7tK+97OLgOrPpY7+DzROVZ44K/hepXewrW4cr3H8f+QeIa
y+6C14lyUp27yg+vZAydg88rlSfq/FY+Zz72ouAzCNX0FBu0Nxdh38p9qiXO
YXW+Kl+6QJqZfsclzh9D8FmG6oeq+yrOvp73uZx88p951+S0Jp+H/y74TET5
Zrfg81blv9uyvaC9fE702YlqyCWCcxnVN1RH1L54hfv8U7AeU363quC8SfVY
6SfpKNWB+2UxS35sSrBmVs3n/wzBHXQ=
         "]], Polygon3DBox[CompressedData["
1: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
3hR+Bd4L/nejN26Ct4bvgacwXqnx3T81uX4UvDfcHkh8vq8v8XdX7Iea7/j8
E7GPg8/TX+AzesnR5mpz+kY1Od3cbk7fpCZHmCv+AG9Vk+PN9WqEzWpioLFw
X/jXNZpB7WCM3KJmTVwbNcEGNWvu2rvma9fsOfeea75qzR4Y0+gT6xT3gHvB
PbR6zR5yL7mH1qzxWX1XHz27xkf1VWPYOTUxzFh2izHaXNiJ7+mDlzbx1tjx
Fjy8Pz6ir1ijXFSjSdWmd5sTa2oWaxc168U1GkWtciT8+5oYZixzD15VE6OM
Vca8K2r2mHvNmPPPmhhoLDQGXlmz59x77smzavace08fvtD54LsuBF9rzK7x
QX3xGnMg/DnnwzhOUCPUaDq1nZpq95rvbH6YoZN3uQ5+CZ6nk/15vZoTnreT
+OtnGo+n6+S7vMoen0is2aEk/pxWM1e1E3/3Ozg3/Z18N9/BeD28k3e7wRoP
nq+T/KGm1D+m7URrnlrzrt1O/FWNp/8N7UT73Qi/qj93kt/UyGrlBTrRty/z
/cbz/bYrid+uwX+wz9bJ2lxtjQPP3kl8uBkeDb/DMTnrP4r7fmZ//QXesdkv
auenmvW3RrZW/hu8X82aPsf5iE7W2j1jfpy+k71kTDY2b1uSH1zj57HP2sna
m9NfhOfuJNf7GX7WjJ3svWvgF+A5OonXzrGx1Rjs3FtjW2vbczi8pmdg7+B0
+C81GlQtquZdSv1dEiuMGep/NbHaWM26ZI0mVhurSZeu0aRqUzXvMjWaQe1g
jli+JgYbi3cwRtZoXLWumnQR+LqSWt6afvYajapWVQMtWqOh1dLG9GVrNJHa
aCtjcm1qYGOla1ijsayVjHFqL2OcsW7ZdrSnNbi1uBp01poegL0ANdqImprf
2l+NOFONprW2VwOqddXA1vZqRrWxNbu1uxp57poeg70Ga4xZavaQucac6t46
v6R2toaeokZD23tQw6qtrUms1dRs1ipqTmtla3y1qJrV2lGNqpb9pYbviyaf
q6Zmt3a3RzNvTQ/BXoIafZ4anzHXmYP1JXtG9o7sGR1bk9PN7eb8leG/ltRO
1lDfdlPjWFtP107tYw1tLW2Nba/LnpO9J3tOx9doLLXWpfAJNT0se1n2uP5U
05OyN6Um+1tNz8velz2ug+CJ9r04zoWPrOmZ2Ts7Bz6iJmY7t1N2EsutGew1
WFNbS5xUUntZgw2u6UHZi1KzLFCj+e01WDNYC6hh1DLWWAvV9LjsTdmjmq+m
5rL2smZbuKbGMDebo609rEHsBVjDWJvYc7A3MbSdXoQ1kb2YeduplfYs0YZq
xPe70ej2Uqwp1e7WcJf2pQa2trMmsbdiz8VaxRxmLpmik9xmjjNXTN5J7jOH
mFt6neSWg+Ez4O+tuXuJAcaC7+zx9JITzA1fGCN6mXPn/mf4PPiv8JnwD8YU
+BD4LPhH1wQ+D74VnroT/awGVguPcr8SP49q9EUL+/nYB/qTW9/2O/ZHP1vL
2gPYfCD61fmzxzMZ9v25/0TOv2z0yYHwyfDX8FnwAfBJ7kX4zF7W3LU/uaRf
VWr04qkl/UFznnvx/aa/YY/R2mxkSe/RnO1e/KDpN9h/Ubs/0+SDk7Bfzvmk
XL+4lz3uXi/whfDR+jzchi+AT4QvgzvwRfAZ8L/ggU70vPtZLTJVJ/WSOdxc
PqyTesDej7noXeNzJ9r+Dc7f41iikxrOWs6cYG7YWR82f5mTsO0JHw1/BJ8C
7wQf7Pxbk8K7wIf4fDVCLznOXPcWfBy8B3wU/CF8cjM/apvR8EnwtvAf4dfh
o+Gt4QPgV9VA8Hbwn+A34GPg38L7w6+oQXrRkGrJl/ULeBv4QPg1+KheamZ7
RZ+V1NLT1/S7rijpH09Z08+8sKSf6BxYf87ZzM3e8LHwJ8agXjSeWm/mTuqx
jexbtLMn1X3nsPfeb2cPLwTPUNMPvrJEfxsTjA23w/vCFV4M/kdJv9gYb6y/
oKRfbY4z111fosf3go/Rt+FT4f6afu8/S/p5xnxj/1kl/Tx7xvb2fGffXY1h
bh9Toj2Mieb+CSWxciQxYxC8W0n9/hbcB/+upL80Ch4M715SL3zWTb3/x5L+
ljWlteWuJfWIMcpYtVdJP83+pfXKgSX1yehu6tl9Svojb8MF/n1JP+/DburT
fUv6TR90U0/vXdK/+6ib/sB+Jf2+1+EWvFNJfT6mm/rmTyX15BfdxIqDSuq1
77vpZx1eUg+N7aY/eGhJv8scYy0wuJPc8wr2CXzE9iX192vwz/COJf0J/dla
oa+TXPZdN/2fw0rqzc+76df9uaSeG99Nf/BYeLtecqa585CS/t+Eburh40r6
+eYMc8cxJf37id3Uh8eX9O/NodYCgzrJrX01/btTSvr1aiJrqUdLtJI5yVrp
rpJcZQ4zl99ZktvUBGq1W0u0ghrM2unuEm2m5rJ2uqdEi6kZ1Hq3lGgJc6ba
7o6SXKq/Wju8VFIL2sO1nh5X0tvVP9Rut5XkUnO62uDbklyv5rCWmaQTLWI+
UguMLdGqala1yvgSLavGUqveXKK91HBqyWtLtJ0a1VrywRLtqqa1dn2oROva
g7Qf8EBJb1JNam36WIlWVcNauz5Som3VVGrfm0q0lprX2vfhEi1sD9Le5o0l
vUk1pNrgmhJt6e9o9tPtydmbUxOqLa4u0YpqVrXuv0q0rJrW2vfeEq2rprcX
80gTz/2NwN6APqfv6UP2bt2T+pY1jvXI0+3UPvrwlY0+0bf1CXsF5kB9xRro
yEZ/WBvpw1c0+snn64P23vUpfdOctVRJjjKXWbMc3ugbaxlz4NIlGsDcqE/a
S3FP66vmXH+bM0abi9Vn9zV6XT1kzPO3DGOssdCcfX+j383lxtQnmvhnrNXn
7KXoU/qiPvFAo6/1FWOEvSd90thhTni+qTfMFeZce2nmAHOxOWG6Ek1jrjAm
T1uieYzV7zf58JcasJccOk2JpjK3GtNHNfWJsX50k3/bjPkrtvdqfhsdz5g/
95ITnmvqF3PFB00+bjH+4Ca/qeUvL4n9L/EOb7aT89VDagL7FfYf1ArW//Yj
rLGtta1h1LpXldQ26vnrGn2pfleDq/2/KdHm1ifPNvWR9dLbTf79kWsH9lKT
29u1R2Ctbo1tL9ma1NpbDThPiSZSG1pPjmn6Oe5tc9Dpzf41N5lz7FWZc8xF
5iR/+zDGm6vUx+c39YX62f6F+t+cZ+5To1ofWA+oXc0R9nLNIeYONercJZpL
7aqGtT9iP0Rta/05sqk39UU1sPWKMchYZE66sKl/zFVq1oVLNJla1pyoP9q/
Nleq8RYp0XRqP3P24039aC5Xs9pbU1OrZc2J/tZlDjFXXt74uzWNWlCN8UpT
f6g91JT2q8xJ5iY1z8tNvedeUXOv1cQLtbg9G3uL9mjs5ahxNmzig9pHjTyi
RJOrnZ9r9Msdfclf9njstduzsPejJrOfZcw39psj7eWb082d1hRzlWhgaw17
QPai1QT2hsyx9rrN+eZeNdPopr43V6jhZy3RzGp7fURfUTOuA7/b6MdxXuvl
8FzNuVZNf8b+ixpeLa/mHV5SM6iF+/oTu5xD51KN2yupGdS+alj7I2p+te3P
NXvHd/Ld7JFOXlJz2Dv9rsb3py7Rw2pk+2v2f9TOnf78PmJMNDb6m7K/XQ/p
5Lfmdxq9/JPv0IsGHyipwdTmpT+x1Jhr7D2EY5KSmLJpzWf62dZwanNjjLHG
GLJxjc/qu2rw1eCxNWN9R9/1+5rYZgwzlk2sWTt9QF8YX7MW7hH3ymEcg0t6
xlti+6lmL+uDv9SLHJOW1AhqfzV9p6RGUOur4YeU1CBq+285TmlnTp3bCfAF
7ewB94J/c+DfWnxa8rcIX9eslWvq2n5TsxauoWt5KMegkhpmc2w/1MRyY7ix
/Meav2UwxhvrrVGGlmhwaxdrkm5JjWCtsr9z1PS3V4APgH/oS82yCtzqT6wz
xhhrRtXYvMd7R9bc+0uPHP6kZi70IX3pzZremD6gL6gBb2/6EWpD7/Fef5Nd
qUZDXt/0D9SWakj/NkSfVVt+VLM2zoFz8VnNXLsmrs1bNXvBd/BdPqyZK/eU
e+vTmrVzzVy7j+HD2llz1/7zmrlxzpy7MTVz55q6tl806+UcOpfWCBc0/SJr
B2PKHU0/S61sf8f+tDWLtYs5z9xnjt2wpobwt1Y1vLXFlzW+oY/pa1/V+LY+
o+9M0p/YbIw31utT+pY5bA3GDupPbDaGG8sHw9e1E7ON3fY87X3eX/K3EOZ4
c71/f+Tvr4c2++enkr9Fsaa4oekvWWts3+T3N0v+VkgN79/iWNOp7a0Z/K1K
DW8t4e+7/n2TfzPl3061+xMfzKHm0v8HkE3ZKw==
         "]], 
        Polygon3DBox[CompressedData["
1:eJwt13mYTmUbAPB3hpnBvMPMa5fJnr0QUdFKKYmK0te+KG20WNtXWiylRUqh
LBWF0KpFiEol2dOeoiRLsvP97uvMH/d1Pffv3M857znnWc5b58p+5/bNTqVS
y0SO2FMhlRqem0rNk7TLSqWy06nUFnYNG8VGijMyqdQu9ih7T36Muix1u9lj
7H3WlqXYP+xGNo5VY7/mJ+frzYazEWJdUSq1ld2Wm5w/rrOeHWLPsc/lp+lb
wfkKClOpWexndjE7gjVn69kBdlDU8Pty2CS2Qt5NXXV12Wwi+4YtEzeqS7Hn
2Rfy09UVqSvNXmbfsrNZNdaG7WWHiZqirr6l2Evay9V1VVdVXQv2F0uLAlGs
rgabr70trqnuGHVZbBxbyjqzDDvofseyz1gnVsh6qGuel0pdwtewmexnzzBP
+yjeQvRU00000u6l7gLxmut2Z43ZhfEc1E/X92zWIPrEc2Cvsq6sPuvBlrFX
2DmsCfsfW8leZ+exZuxitprNYOeypuwitoq9wc5i9dh57Gs2le11b6Pln7i3
+eIav28fe5ItkJ+groy6A+xptoidxMqx/ewptpCdyMqyQteYw35jl7HGrIjN
ZRvY5awJq8rmsc2sD2vFKrDZ7Fd2KWvEyrLX2Dp2PqvNyrBX2VrWkx3OyrM3
2S/sEtaQlWPT2HfsAlaH5bPpbD3rxeqyauwj9g/bIm71DNJsJvtJ/qO4OZPU
fZib5OGlWRX2PvtLfq3ztXS+DHub/cGuYkeyiuwdtpFdzZqzXDY53iPrzmqw
fz3T+9nLbFJcx3yrpO5dtkneO8aXusrsPfZnvDPWouRdPsHeZLPF70XJfczI
TeZVzK/drAuraxycy7/Sd0rJHJyQm4yBGAtb1BUUmMfqTo93rG5arBviO/l+
NUPYFfL2+pZV14Q3FY09l8V8iPZYdVXU/WNuLGOPxBxgR7BD7Gv2MHuVNWAH
2VdsGHuF1WcH2OnsNFHftb5wbLBrfCt/XHuWuiPV5cmXi1FsJmvOcuWfisFs
DKvMtjjfO+xq9ijLYt+xt9iVbCg7ZAlay+ayK9iDsYaxNextdhV7OBZmfdex
91gfNjrGBfuBvcuuZY+zUux7tpANZE+yimwz+4aNZDNYM5YT81x0F43c7zeO
3eV+v5QP1Z6qrp66/fqW8o5Snn0d/gJ7qOR9lGPHsWPFS/p2YPnsePlb6kap
y9O3DGvMprCRMSZZHmvEJrMRse6yLFaXvciGshOcL83as7fZ4+wHsTjWWNc7
Sftk8X08f5bP2mt3EEvE7ex591FV363u40TnK3C+Dvwd9oSa2eIy+b3q9rvN
Veo+Y3ez8aymup3sC3Yvm8gOZ7vY5+weNoEVs//YUnZfPA9Wi+1mc9jl7H52
wDVWs/fZdewplqPuRzaPXc+eYbnsp/xkf/u+ZD+JfWWYe1yl7lntd9W1UVde
fqS631ilGAvsenYU2yCvLJ5g18W8Em/E2qzvWayBfIV4Jp4xO5ql5SvFmHhO
rDUriHXXO/pA+2Te0nMskq/mL8o/UtfescryZq67kuWI+9hFbK2YIl+s7lRW
M/aU2ENEY/XLHbvbva2TT9Veoq6jumJ5z9gLRRN13zp2j7o18vHaH6vroK6K
/AJxvmiqbkW80/gW8Btz/daG8knqhpeMjxgntQuTsdRf3Ykx3kStwmQs3cbm
y2+OZ+4a5fX9w/v4hN3CHmEV2EZ2BussGui71LEh+i6Q36r9mLpCdZvUfcBu
YM/mJHt47OUfs37sIVbAfmcfsb7sAZZmG9giNog9zSqxv1m32NNFQ9dd5tid
rlusvVD7X3X91B3r+OFsEdsZeww7rmTexPyp49gSxwboe6r8FFGXfcYGsq6x
p4sj2NfsDrZe/lZ8Bzjfec7XTN5JdBT11H3u2CB1nbSr5SXrQ6wTrVi293FI
u7YYp++D8VvUZdSdwj5kz7BdYqO8bKyVrK98t9gkLxfvhPWT7xF/xtyPtZLd
nE6+1/aV7BOxXzzvusew/drNYm1Td6+6tuyAvHnsA+y+2BvFr/LSMX5YH/lO
8bs8L9ZodqP8L7FavtszGMAuiX7p5P0Xud6Z2l3EL+lkLynMJGMkxsqmdDI+
q7Fe2hfGPYhVbJfz9Xe+i+V7Yy+MvmI46y//T/whLxNjht0k3xHXkZeKucGu
lf8dc06+1/kGscvkm8UatocNZJemk2/Hhp59Z36G6OA37bDH36k9XV0jdduM
tf/YQ2wua8H2sp3sQTaHHcX25CffFg/kJt8G8Y1waSb5vu+fm3xrxDdHD7aN
DYj1X15b3z/13c7uYNNYQ7adVff7PmZb2Q2sjd+8Raxj+9hgdrl8u/iZZcc6
wq6JPJ2MxQqZZF+P/X2b+IllxX6nrrf8MNf4RL7d+W5ibdmZrI7ncg7/kk1m
W8WP8lTMNXZ1nCudzJXymWT8xzxop2+2vq34G+oejn7pZJ4VZJI5FnMt5kf1
vOS5x/OfGnOQVWIdY23Td2zMK1aZdWLz2XMsbR7VYiexeWwSO15dDmvDZrFH
2WmsBjuTLWIvsv3xHOQZMZINjHGprlhd11h72USWco3d8uqxJ7G7Y31QV1Hd
qbHms2fZwRiX8qqxn7E75cepK62uNZvJHonnzw5jXdinbHysm6wmOyv2BjaB
ZbnuHnkNMYbdww7FmhBzJtZAdlesX/oW5SX7UexLY9ixrBQ7ms1gw1gZ56sQ
/5diTLOnYj6nk2+eyplkD4q9qKa+C9gO46CvunYx19LJ3lQlk+xBsRdtTCf7
UNVMst/EvvNbOlkfK2aSNTnW5tbOt5M1EBOd73a2IZ2s1ZUyybdXfIO1VPc3
K469Qd0t6eS/5Obc5P9m/O8crb4V26FdT4xXN0RdWfdW6N5astfZ0yyfVWXt
Yr6yF9jR+v4rry8msMHsQMyb+M1iFBsk3xfjPNYwMYINkJdzvirO1zbmOhvH
8pzvFfmamM+sFvs/uV8bLg==
         "]]}]}, 
     {EdgeForm[None], GraphicsGroup3DBox[{Polygon3DBox[CompressedData["
1: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
iPhdwXn+E3h8sJ4vEO3/qndOot8NqMcpFt2vymMvEd3fyoc/M7ov1XuHM6J7
V70jKBKt6fXuII32rNWXlYzuV+Xtnx7tX6v3idEevXqrs6L1tN5xFI/umeXz
nxOtufUe5JTodxjyPc6O1tx6VzIqWvfLlz4p+v2B+qzB0b2KPMY1qd9d6Z1U
reh+VX5s8+h7KN+vW7Rekbd5Q7RekQ/ZMLoGy0fqG13n5Gf2ia6L8lRPje5L
1at2j9Y38lSvis698r0HRNc5edrto2u5vMr+0TVbnnbPaP0kb/b66For3/Xq
6Lwtj/Sa6Dolz7NrdA2WJ9wrWrfJ7/0udU+lvuam6Lou/7x6tF8jz/mKaB9H
7+8aRdcO+eH9omu/fPXG0TlE3vjO1D6dvLum0XNJ3mmb6Boqv/eyaK9EWrpJ
9DyUl94jWi/KT24RPffk161I/d5O7w2vi9Zb8rFvjNYHenfQLHqNyGceEd3b
y9vsEK2f5OveGt3/yzPvGK0/5BsPj+7/5cPXjfY49F6gbbR+ktddO9qzkG/f
KrrO6n3E5dFeiXR+jWj/S7596+j1q/cC9aL9Eb1fuDZaY8kD7xyt/+Q/D4zW
MXrfcUu0jtF7iirRvozeDy5L/R5U7z3bResz+duDovWT3onUjPaz9K7h59Se
nXy8+tEaXe8m/g+E0wiD
         "]], Polygon3DBox[CompressedData["
1:eJwtmnfgVmMfxn8e43xf3vs4+qFIEpVSySyjQoqWUIr2TnvStveWUq+9s0Kl
SINCkUQl2Xvvvdd7fVzPH6fO59z3c57n3Oe+v9/re92/Wv1HdxpVqqioGKh/
ttL/RWVFxRZRUfGE+KtUUdFVRxOdn6xrLdS2k/hl8SYdHfOKiu/F7+j8dR2d
xTuoT0l9V4i/VluleCvxk+Jvxe3VZ7nO/84qKiaprZ14mfgv8UTxW+rzoPhF
fWayeDe19xDfpPajxVV0bKm2lbr2jfrOr1JRcZLaKsUHipPa/xQvFX8mrq2+
h4i3EtcW/6H+74tv5TvEVXX/Vjq/QNea6bMPqL2JzrfUtT3VvpM+/5GOH7eo
qNhRvEHtV6v9eLWPFu+szzfX+Vm6drg+30R8hfgFcV/GS9xMfKb4MHEj8Vni
FeKu4h3Fh4vPEB8qnqv7t9L5H/q+3XT/fdV+ttpX6trJar9N18aJZ+g3dRMf
K95b3IbnFDdW/3PET6r/KeJq4tbiC8XNxSfpeFzn03TtRn12B7UfpPPTde1g
tf1Hx6e616O69pHat1d7I55V7QeqrRDvKx4jPkh8kPp20fmOunaQ+j+k399G
vJ24Ae9SfIx4W/E+4lyfb6jzUbp2gD6fxA3EI8X7i+epf2edV9G1/dV/gbiH
eBfxoeK/dRQ676RrZ6h/Ze75OVHcVHyI2qvRV79rP/F2aq8rPlXt+4q3FdcR
DxY3Ev839/gNETeu9Dvnu/cKz4XTxCcxPuJWas8rPVeZs5+rbSt9vrrOe+na
3mrbRlxD3EdcXxzimuJ+4gbiTLy7uK94H/FafcclOj+W36j7ba323XTeW9fq
qf1htffU+a66dpjanxdfKm7DbxYvrOLv5jccLl4mHiveT3y8eLl4HGMrPjH5
HTCWDcLvhjGuJT4gPPaP6hio/nvp2tHJ75Dvbhh+t6yxPVlr4bV3ZPJ3HSVe
JH5Ex4DMfVqqrb6O//Du1V630r9xH3Gz8G/nmerz28PP+gRjzlwjnqitQuOx
k867ck/mb3iuMmdPUntJ7VV1foqu1SZ2iXcmXvEMzGdxbfEgccNKz6F6jGV4
bj2uYzxzkTmlti1zz59uulZH/QeGYwMxYqDa/0me68z5WpV+Rsa2ZfjZGTPu
dXR4LFuFv5vf0Cp5jPpnHnPGbkfdY2v1eUr8nXhpFa8F1kR78aDk9TtYfVqq
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
         "]], 
        Polygon3DBox[CompressedData["
1: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
7Pq/ALqUKGM=
         "]]}]}, 
     {EdgeForm[None], GraphicsGroup3DBox[{Polygon3DBox[CompressedData["
1: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
iomtQufQWkedQtdstEY6h645aU11gS+OXNu5Gp4Wue73n9A5h9ZRt9B6WXpv
buj8Q2uwq+ZNxu9E8WN65Pre46HzD62jnspfI/vMysprI/uNSqHzNs3hispB
Y+fY8s/KOfph76c8O3aeL58/LbaGkP5U3qM12zd0zqd5W0E5buz8X7HgROx8
WPErK3Ger7h2PnYdQbGjc+Qa4FjluLFze8U4xRvFnVKhc1Zp+JLKm2Pn6qqd
al4prjUKrR0UB5sq/46d58sHlo5cA2kSOj+TrusNz4tcn3w+zRekaZ8LraMV
a2qHzv8Uo/uE1lu6VmO9k9Axuqz8QGzfpPVbNHL8rYv9n9i1Ca33s7FrEPI5
p2PXF+STk8gxuo5iVOjYWob+TOy6g/yDcnFp3RL0m2nHWHd/ZTuvlQ65Fluv
yDXASfC9kfXzM/CNkeuNI6T3Imvsp+C7ItevHgydJ2g9PqxcKbJufxnuHbke
OzldF/KN40Pn5dKKDdK1rPf2bOj8TJrtodBzUXPyVrhH5NrXBLhn5JrYRHhq
5NzhldA5gWpis+EBkWtosxRPI+t56Xr5LdW3x8G3Ra5RT4fHRdb2T8OdItdO
75afj1zHW6i1HzmveRGeHDlPeQm+JXItdzTcP3KdcCZ8e+S694zQuan058r0
Xck/zIO7RK693wMPjVwnXAAPiVzTmx86/5bvnRI6h5a/nRra30sjlQud90uL
XpfOMenhNqF9ut75faFzQdXYl8NjItcbl4WueUiXtg6dfyu3XRy6nqE43j50
7UHx98bQ9Q9p2utD1/MU92+GR0eu5y8NnbsrL14SulaheN0hdF1ccfam0LUc
aYmOoXN65c6L4Jsi17RHau5Grrf3D53HK7ddG9qH6xnbhs4JtE/RLnRdQXnx
o6HrCspV14SuByhHfix9J4p3Q0LHbdXn74DbRa6xD9W6irxHcLvWW+T6/0Ct
w8g1/wHyC5H3DgbBN0Su2w+D20eu5w+Xv4u8jzAYzpW4xidN2zFyrX4M/y+O
rVGUd/8YWyOqVnAwtr5UfeBAbL2oXPvCxLUL7V/8FFs7qm5wNHYdUBp4c+x5
r3xTtTfl3bH6xL5MezplEtd8tDdRInGNRXsZRRL7IO2V7I+tX5X//hy7VqJ6
/i+xa4iq5++NXWdXnr4vdv1dufw3sX2B6hJbY68Z1Su2xV63qlE0SrwOVQO/
KvEc0r7Sb7H1gfTzr7E1hPYUfo/tv6RdiyeuZWkv5o/YtU7p3sOxa53KCxom
9kGqyTdOvCa136Q6gPIU5SuqXSj/6gi3S6wVlOe2Tezflbe2SLw+tb/WMrFf
0N7ZZUm6/uHCiX2x9lauTaxjtG/VKrGv1B7c5Un67eHrEusG7Ze1Thx3pfGu
Txy/pf22x9aLqm80T+wXtK9XO7Ff1j5Xg8S+W3sNpRPX7rRv1SSxP9KeQtPE
PlH7a/kTx2/ldHWS1AfB9RL7We0j1E/sZ7VnsSTN3y9XTE9S/wiXSlw/1N7Z
lUkag+GSieuE2i/bEttvql50dWK/oz075SrSctXpg8R1XuUgFRLXD7W3pVqy
ahqX81uVxGtP+1lVE69J7fdVS7wOtXeZJ3HdXzndBYn3CZRDVUysObVHptqz
6iERfc3E61Z7K8pdlcPmp784cd1S+2WVEvsC7dNpH0D1iisU3xNrU+2JFEys
EZULF0isk5Q7b4wdL1Vb+Dp2jFH94XjsGr1y5DCxDlBOuil23FW9QjV4afVq
9NmJa8rK+yon9ona48uXON4rB7w0sd/Uvqpq6qopXcZveRPvnSh/rJHYx2n/
qHzi+rD2HGsl9lna+/s+dg1UNa7dsffYVMcrl7jGrn3JsonrvdqL3BV77031
ve9i149Uf9sRu06h+lvRxNpUe8faT1B+UZV+T+z9PNX6jsTez1COr30G5YOX
6B0k9qHajz4We/9D9YScxPpe9Yr/BdJGFHE=
         "]], 
        Polygon3DBox[CompressedData["
1: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
1agBk5v7NtHyZ/k73VpMP15Yz6Z3a+FLevHIeuXx4fN6ifHG+uH9eAX1W28x
QT9eTw+tlx7sx0vrKfQWl5Ssx729+Kfp+vG2eii91LT9eGM9r95Xv6kf0l/p
L79vE9uMgcZCawq3ValBWWtYFX4YfhiepYnnd73mKskFjKHGwtlLYqsexfvN
UeJdjHGu95wlsc8cQu8xT0lucTV/q3L/VTl3K2NXDSZ/Mv8z37sAXhxeXI2E
T4fnhOeEz6mTo/i85i7JXdRca59LlTzvhfl/I+Zv7G8ajGboDSYq0RJzFmNT
W5LLGGPMTfslsUcN0ntMXKJNapSxZ7IS7VJTjIWTlmiNv2FjeLaS3+Zv0MvN
XPLb1CS9ziQlWnU/fwfy3Q7k3BjG/gUfDx8PD7DeT/idqrwjDfwPeCd4J/id
OvvDGD2687OH+A6zX74eEv95vDHJ2Mj8Y+vUYA+okvNbm1WDN+r0S23Wo63S
+SW9mzXHq6p4FGuRU8I7VtGgx+r8qXuus2NH+s6qHYzvX+f9NGbpmfXOhzH+
S+ePd/D3tslFzEmsRT4Al15qTNaafGY+O9f0ceYvUnLO4yf4e7BN7mYOZ+3k
/jb31oNbG/l3m9qaNT1rew+3+W7WJKxNPNomlzOns1ZpTvMZa/dklVznvjax
xxhk7fCzNrmxfsf393N47V5yaHPpMW1ybXNsc239kv5ZP63/s0ai/9RP65+s
EVsrtkb1a50albUqa8q/1KkhWUvSr35dp2Zm7cwa0zd1amzW2qypfVvHw+pl
zQF+18Sz6d30g/pda9TWqvW3P9ep2el39bP6ST3Ezt375fq+0WZ/WLOxdvNq
G2+lp/J9f73N3lTj1fr/tvFaejr3xyttvJyezPf9hTZarqfR27wEL9CLp3M/
vtzG2+kB3W8vtvE+egq9hR5Zr6xnH79JjmOsNGbq180JzA30rEOb5DTmNuYM
YzfpIdhL0D9/WqcGaC3QHsPoOj0AewHW7N+v08PQb+uP7S/Y87D3Yc/ggzo9
FHsp9kQ+rFNz1B/r99/t/P5B3f52v9qzsHdhjvRWnRzJXMkc5k3zxTa1FD2P
8UcPqBfUI71Up6Zmbc0a5md1ch5zH3OmN+rUMK1lWoP7vMt31Fo11/xNz6X3
MkYaK/Wcek892Wt1PJheTA/2Qh1PpbfS4z1Tx+Pp9fScT9XxaHo1PdhzdXI2
czc97evwc21itTUi4+sTbWo/1oCs9T/eppZvjd9a/9B+1to1N9cZBp/WS45j
rlP1k+uYo5ir9PrJzcxhzGXG6metXXNzryHwqF5yOHM58x+9nJ7O/Gjcfvam
e9Rc1xqktUjrMcOtJ/WzN9wj5oZvt8l99dPm1qWfXMmc0NzwtzZr6xrbm/qy
Ta3Nmpy1uS/a1NqsyVmbG92mlmbNTf36ro2emnOYe3zbJvaZc5h7fNNGT8w5
zD1+hffspedn7+/rNrUYawjWEr5qUwu0RmOtZpx+tMUc31xfj+7vGSzZXx+1
0WJzJPX8kza1B3MW9futNrm+NQBrAY+10RZ7MvZmzMHMxcx57u/itXqlZqgP
b7bRV2sG1g6ebtMLsmdk7+hJeMZeejr2dp5q09uxZ2Tv6Jk2tV1rvPaOzKnM
rcy37+XzPm6T21kz0X+YM5o7miM+WCdnMncyR3oUfreNH7DGYK3hwzb6rgdR
rz9tk6tZ89TvvAMv20uNwVrDe238gTUMaxkftPED5ijmKuZg5mLmhI/Uqdno
RfWk98Hvt/ED5lTmVuaM5o7mvA/Xyflc7wlK4uUjzJ+0l56RvSNrEu6HzUve
b3tIrtf33bg9LNf3h07/rMG7377s9M8auev5Y6d/+umqFw+uFz+xibfUY+q9
9RTq6/NVvIaeUK94VolXtKbr/v6q03M1SC3bu0SbrOGqD6Xk+Y3o6lVHlnjt
XtdP2LlEu8br+kt7lcT6sbv+0a4lsX1o12/ZpST223NUn6qS98keo/ozrOT9
OamJlzVHMFfQQ6pfTYm3HN7Vlw4piQXWTMzXDijJT9QA49OaJfvbGGIsObAk
nxns6iUHl+Q+1jytBb9epRZqDcd8bv8Sv+8aX1SSw7j29uyNZ2uU7Gc9oft1
gxKv6HqZixxRkh/oKYxfG5bEOzVRL7FPiVZO3tW3jirJVdR449mWJfHVmo/x
a4uS+KlmGdtXKdEya/jmgoeW5CvD+B0v9vJ9LyzJscx9zixZTz2d/uqdKl5P
j6G2bl3iPcwxzT2PKck99Sxq52YlXsaatH7w3Sr+Tc+jFm9S4oX0PGrvpiVe
SA9lbr9xibfy96ntW5V4GzVQrV2vRBvVdLVqtRKtV8PVuj+VaLvPVy1cq0T7
1Vi1ed0S7VXj1eLVS7RfD6HWrV3iLfQQavE6Jd7iH/yOl3rJCa0ZqPnGnhVK
vIAaaWxatkQ71WRj2dIlWq1nVz8HSry8PXD905ASvbQHrl8YqyT+mhOrXSeV
7Fd74vqpcUviqR5Tb7BtiffUY+oltinxnvZM9Fu9En3Wwxvfx1Tx9noMY9Ny
Jd7DHrD++qPOr+shjP3LlHgLeyT66487f26PRb/+ahU/vXiT3oA1fWv7eiBj
5fIl3siaufWVUSX5rjUcvciqJbUdewLfsb9ertIrMEdUv6YsyR3tsdkbGF2l
92YNSn2YqsSv2pN7ED69pH5jT1K9m7rE79pz8ft9WqUXY81bvZ62RP/0dPqN
lUpyM2vs5gfTlfgDa+LqxzQlemf8NDd7oUqvwJ7h9/ArVXqJ1lzN93cv8TfW
sPULI0v0xpqX+j5Tib7Y47E38kGV3o89IHuRH1bpDVkzNx/5fYm+2cO0V/V5
ld6m+83ex3tVekP2eKy9/K9K78eepL2cz6r0Ku2pPACfVlIvs4ZgLeGEkvqY
PY374VNL6mf2UL7k9z1bpbdiD+yfjJ9SUt+yhm49Y48S/2SPyt7Q+1V6V8Z/
vf2eJf7Lno+9p9e698Wem72at6r04syZzZW3L8mlzZmtR+1Xkkvbk7KX82aV
XpU9NOt9h5XUZ+wh2Qt7o0pvyZ6l9cHDS+ovemi99I4lftGe2RfMf7pKL82e
3+fwU1V6gfY09H8zlPgDa27/gs8uqcXZc/qG+S9W6UXZE7Rec1xJvc/4ai3r
jJL6oD2zb5n/UpVemjVY/d2MXW3WHqX1n+NL6ov2cKyv7VDi9+wRWI/arcRP
2+OzXnlsSb3u/xGU7G8=
         "]], Polygon3DBox[CompressedData["
1:eJwt13fUVMUZB+BLE9jIt+DSlS5gBRs1UZpUoxQjTenSlSpVpUq1UIw0C6DH
mkQFo4AxIhClxd4VVFSwRHpHSp73XP54z9l55jezd2+ZuVul99AOQ/InSfKB
KqRaZpPkp6JJUlgVKZAkQ/OSpAXbpV0onA1hZ4onyc2ZJJnKprJfWLZEktzN
/sZeYCWMa6q+1T5aJEnysQFyH6uP1HFzNJavWjJJWsn9LFc0Sm6Y/n+rN9Q+
uTpyZeVO+XyDz3fLjVeXsvdl3lNH9V2nrzK73nw/6s8Xx22+O/SXcnzT9T/L
nlOf5ZJkA1+vDhhbT195YyvJPeLzazKr1A65jTLvqENyDfRdILdW+021l13D
yrC3tf+jDrL67Hz2ofYH6hhrxKqw0z538HmK+ac4vl/1N3fMO7ULqHPYnaw1
+027mMpjo9g69Zbab4665ihnvmZy38kcd57zyw3UX8Dv6K1/Pp/PjrFr5b6M
+SPrwvdijdg27f3GnmZ9WXt2klVQFY2dxt5V/1VHfO+15q3ke8/1HWMz6fmM
89qCFWPj2PNxno3NmutmdVq7kqrMZphns9qkDpuvoXwFY2+Q2yOTVcXlxujP
M994/S+w51lxmYR1ZbPZLFXX2Bv5IZ9LRcndY+xN7LB26Sh2L8tvbC9j57F5
7DgrxAayxWyRamS+a4x9z+d9zste56WTXBG5YXLL+DJji8h0UKfiPMXvY9Pl
2rIj2mWi2ATWjh3VLhvFJrLC5htsvkfZEpZP5hw2iC1hi1nCOoazmqwme8jY
onLD2XK2nBWWuUUVYjVYDfag3KfqE3XCeW6ir5rf1kDuU5nTftthv62b/qvZ
u2wP28M6sobsM5aoI6w7+1x9pn43X1PzXWi+unIfyBwydr9cF/0n9TfTP5KP
dCxfs/pyH2kfkzsodyurxz6M+dkB1pU1Yd+wwywxtj/7I/uc5Y/1RK4H+xP7
QrtgzMl6sozzMsL3PsmeNLaozB/YSPYUe4pl2F9UQVadVWcPGNuYbdc+6HvP
mK8fy2fsrXL389ly+1kbud3xTMa9ykazP7O9cd+qEmwsK2hsH2MfZg+zE6yO
3PvaB3zHPt/ROdYW5+ogq6aqyo1n98r1NHYym8yOst1yV7HerBdbw7pl0/Xn
pVy6Xsa62T2brpkv59K1ItaM2ew+9rSxzxh7QTa97wdk0mconqV35O/n09gz
8UzLVYjjyKbr48pcur7FOtcrm65xr+TStTHWyF8dX23Ww9juxr7GvmRlWSvW
ii1jX7FyrDVrzZazB833V/Yqe41dqt0jm65xK3LpOh3r9c/GXs66yd0Wz7X7
r182XVvX5NK9JfaYr+XKszYybcz3JHtA7kH2InuRXajdP5vuQ6/n0ucknpcN
xv4qUz5KbnjsE+wX7XJRbESs/+wT7VOu5Um1O5fuUW3MN5jfod6MY2J5rLb2
FbEes1m+dwJbqr3UfKWy6Xo6JpOuQbEWvSv3P2OvYD21e8it8r19s+letzqX
7i2xx9yeTdflVbl0P4x9cZ2xPxiXU+cZOyj2DvZ9tOM+jedGvk82XYNf9XmL
zFbVO5vudf/Mpet0rNddWSlWN9Zc8z3OtpmvIruJ3ahKuh5fsDKspXZLuSfy
0nOQZY3YdfGcxN7NirD62vXl5sS1Y8VZ41iH2UL2HavKOrD27Hm2g1VjN7MO
7IU4RnZC+5IoNiX2GHZc++IoNinOF8sY25A1ZPPjnmYH4vlTVdg4toUlmfSa
xbXLOOb3WGFWT7ue3Ny499kx7Yui2OTY81lRuQasAZvHvmVVWHvWjj0X+yor
wK5mV7HZ7HtWnd3CbmEvsh9ZTdaJdWIvxbrOSrAmrDFbFPcvO481ZU3Y4rPP
zGWZ9PrEdXrONf2clWYtYt+WWxr7O8vPrmJXslmx1rMca8aasiVsO6vE2sZ8
7Jl472H749zFHh/rtHO1leWTuzLue7mZct+wyqwda8ueZT+xS1hX1pWtZDvZ
Rawz68xeZrvYxawL68JWsB9YDdaRdWT/YF3cpyVZHVaHPcZuZRewa+P+Y0+z
tcZ+Fdc1ivVh98S6IzeJTWJHzt7jFTLpfRv372PO3zC5G9kQ7SFyP8oNZzex
oWwo28nuZt3ZRDaRHT77DO4omq4ZsXZMMt80ubvkHuEL5M7NpuvLx9onZH5X
E+QqWyMWym2M9Vrusmz6zrqAvcPejmvi3M/k97IntJ+QK6k9Q41mC9kilqc9
XY1iC9hCVkx7qGrNBrHBbIdjruA75rP1bD27SGaIasUGskHsO7k7WAvWnw1g
37Df/I4rWS/Wk61mY+K9ho1hY9keNir2ZTaKjWa/sV+MrcW6s27sVTY63nXY
aDaG7WZjWSc2lo1je9l41o1NYBPYITaOdWbj2Hi2jw1i17O+rB/bxgaz5qwf
68+2sxGsLRvGhrFdbCRrx4az4ewndhdrz0awEexndidryQawgexbVtE5fZht
YBvYxXGPsnlsHVvHarIqbBHbxDayy9kUdSd7iD3ECsb+HteOzWVz2Tnak2OP
YzPZLHbG905kfdl0NoOdYpNYPzaDzWSnz65/Z84+z/Fcz86lz/m+s3tH7CH3
sQHZ9N3yX7n03TDeETuzYqxWrJ3mW8A6xR7HLme12CN56X+Q+C/yeNyr8T5W
In2/75JJ/+PFf7115r3N2AtZc+3mxv7duNJyM2MfZKtYRZkybFa8C7DVrBIr
y2bHHs/WsMqsHLs/jpm9zqqw8mwOW8vWsursfDaXvcXeYjXiONijbCvbwq5g
VdlitpltYrVYNbaEbWGbWW02MJu+N7/hN30R70dqj3N6Nesj11vudZYzdkq8
c7GVrJxx57HJ8T7EVrCy8Q6hbmfT2HR2Mi/9rzsjk/5/jf+xN1gPpsbzyuZo
z5ErpF1CblK8J7KXWWlWnE2M90n2EisTaymbGu967BVWns3Npv2RW6pvnvb/
AanWIK0=
         "]]}]}, {}, {}, {}, {}, {}, {}, {}}, {
     {GrayLevel[0], Line3DBox[CompressedData["
1:eJwt0U1LVFEcwOGTlemMo61FlJYithbRtUitXIjkLnChC1u0sUWCK/sE+QlG
fNds3Iu6N6d8LacR0XJ8RbNSRHwOCPd3nwv3crjnf569ftPW+yCEMKC/is/n
2lU6EcJTfmYNF/icK2xmni95yk7esoepZAjvWMUPrOMQGznMVs6xg0vM8hf/
66v+6Egjj0L4xlGucoxrHOc6J7jBSW5yiluc5nfO8Ac/cZuzzMX98CczzPP9
4xAOeKUdXepYb/1Xll1cZDszbGGaDfzIWg6ykn1Msps35vEqrsUXzLGJy6zn
PKv5xfd7/KcKXegwzvNhCOVxHaZ4zbI4JybjeTER58XSeG4s4Rmf8ITFcZa0
xVCgLYffdIV9tyL2e1m4P/87AMVRVw==
       "]]}, 
     {GrayLevel[0], Line3DBox[CompressedData["
1:eJwl0jtPVFEUgNE9wMxoACE22Bhigw2hYmhMoMFeeQgqNNAYC4w0xG4kQAjB
ODEgxgIaWiio9ReA8pCHEGKMozY+QEtRxnVD8WXtW91zds6VgQcdQ6mIGNeR
4Q+/8VA7lRHnWeQl/uZVpqoiWljL66xnF5s4yFa2adh8l495nwU+4gInucw5
vubF6oiP/KU1/dAX3chEvOFNvmUH19nJDXZxk93c4i2+Yw97tW2+zR3e4W5y
Ju6xj+/Zz302ZyM+JffUgX7qq1accZWLfMVZLnGC8xzhU95jPvkfe/XQfI0D
bGQnL7OdF5hjyU4beMw61tjDh+Rb5/Rdn/UyHZHlC2b4nGnOsILPWM4Cy/iE
00qZpxjJzlkyTvCUY/zHUf5lnifM2UMxdfYe/gPy6UqI
       "]]}, 
     {GrayLevel[0], Line3DBox[CompressedData["
1:eJwl0rkvZlEYwOHja2xjpx6iH6IREURIZGLXzohEaYQaCZWEvwDRI3pkevu+
ZCJihsY6VPad54vil+c9uc2977k5rZ1NHTEhhH6dRUK4djjiH2UlhRBhLtNY
wGyWM48NLGULa9nBn+xln36ZR9nNSQ7yN4e5yHHucppNySHs8b/mdKIDxcaF
MM84LjCei0zgEhO5zC9cYRJXmcwUrZlTuc40bjCdm8zgFjO5zQvuR9WOTnWo
K+84yyNORXfDMc5ziDMc4AS7OBL9RrWZe/iD7axhM0tYx28s41fmM5X19rDL
c8XoWP/UGOuCos/57o5q+cZqvvI7X1jFZ1ayQk/mcj6yjA8s4T2Lecci3rKQ
N7y0h7+Rz//hA8AMTBo=
       "]]}, 
     {GrayLevel[0], Line3DBox[CompressedData["
1:eJwl0skqxVEcwPFj6OZasBKWSrG24Al4APECWCMeQFfZmva6YiseQObhDiTJ
mMwWknm4uLfkc7P49vmdxf90zulf09HT2l0QQhjSVVEID4UhnHFH7WUh5Kw7
WWrdx2rGWM9RNjHOFs6ynYtcUpf5lP285yB/OMaS8hAmWck5jnOPN/nvdaGj
/By1F5e4zGWucIWrXOUa17jOdW5wg5v5zAkmmGSSKaaYZppb3OJUaQj7vNW2
LnWsBmdbYC1nWcE4IxzllzvFeMc+nvBYneZ5tnGGzZxgI4dZxwFWsZdRjthz
l9fKevtzHubfyxl/rCP8ZjG/WMgMAz/5WxLCB3PM6t38zTdm+MoPvvCNz3zh
E5/4yGnvcFD0/z/8AT6nXZw=
       "]]}, {
      Line3DBox[{2725, 3023, 1205, 2724, 4390, 2932, 2726, 4391, 2933, 2727, 
       4392, 2934, 2728, 4393, 2935, 2729, 4394, 2936, 2730, 4395, 2937, 2731,
        3988, 4480, 2732, 4396, 2938, 2733, 4397, 2939, 2734, 4398, 2940, 
       2735, 4399, 2941, 2736, 4400, 2942, 2737, 1401, 2943, 3024}], 
      Line3DBox[{2739, 3989, 4481, 2738, 1219, 2740, 4401, 2944, 2741, 4402, 
       2945, 2742, 4403, 2946, 2743, 4404, 2947, 2744, 4405, 2948, 2745, 3990,
        4482, 2746, 3991, 4483, 2747, 4406, 2949, 2748, 4407, 2950, 2749, 
       4408, 2951, 2750, 4409, 2952, 2751, 4410, 2953, 2752}], 
      Line3DBox[{2754, 3992, 4484, 2753, 3993, 4485, 2755, 1234, 2756, 4411, 
       2954, 2757, 4412, 2955, 2758, 4413, 2956, 2759, 4414, 2957, 2760, 3994,
        4486, 2761, 3995, 4487, 2762, 3996, 4488, 2763, 4415, 2958, 2764, 
       4416, 2959, 2765, 4417, 2960, 2766, 4418, 2961, 2767}], 
      Line3DBox[{2769, 3997, 4489, 2768, 3998, 4490, 2770, 3999, 4491, 2771, 
       1249, 2772, 4419, 2962, 2773, 4420, 2963, 2774, 4421, 2964, 2775, 4000,
        4492, 2776, 4001, 4493, 2777, 4002, 4494, 2778, 4003, 4495, 2779, 
       4422, 2965, 2780, 4423, 2966, 2781, 4424, 2967, 2782}], 
      Line3DBox[{2784, 4004, 4496, 2783, 4005, 4497, 2785, 4006, 4498, 2786, 
       4007, 4499, 2787, 1264, 2788, 4425, 2968, 2789, 4426, 2969, 2790, 4008,
        4500, 2791, 4009, 4501, 2792, 4010, 4502, 2793, 4011, 4503, 2794, 
       4012, 4504, 2795, 4427, 2970, 2796, 4428, 2971, 2797}], 
      Line3DBox[{2799, 4013, 4505, 2798, 4014, 4506, 2800, 4015, 4507, 2801, 
       4016, 4508, 2802, 4017, 4509, 2803, 1279, 2804, 4429, 2972, 2805, 4018,
        4510, 2806, 4019, 4511, 2807, 4020, 4512, 2808, 4021, 4513, 2809, 
       4022, 4514, 2810, 4023, 4515, 2811, 4430, 2973, 2812}], 
      Line3DBox[{106, 1094, 107, 1095, 108, 1096, 109, 1097, 110, 1098, 111, 
       1099, 112, 1100, 113, 1101, 114, 1102, 115, 1103, 116, 1104, 117, 1105,
        118, 1106, 119, 1107, 120}], 
      Line3DBox[{2814, 4024, 4516, 2813, 4025, 4517, 2815, 4026, 4518, 2816, 
       4027, 4519, 2817, 4028, 4520, 2818, 4029, 4521, 2819, 1294, 2820, 4030,
        4522, 2821, 4031, 4523, 2822, 4032, 4524, 2823, 4033, 4525, 2824, 
       4034, 4526, 2825, 4035, 4527, 2826, 4036, 4528, 2827}], 
      Line3DBox[{2829, 4431, 2974, 2828, 4432, 2975, 2830, 4433, 2976, 2831, 
       4434, 2977, 2832, 4435, 2978, 2833, 4436, 2979, 2834, 4437, 2980, 2835,
        1309, 2836, 4438, 2981, 2837, 4439, 2982, 2838, 4440, 2983, 2839, 
       4441, 2984, 2840, 4442, 2985, 2841, 4443, 2986, 2842}], 
      Line3DBox[{2844, 4037, 4529, 2843, 4444, 2987, 2845, 4445, 2988, 2846, 
       4446, 2989, 2847, 4447, 2990, 2848, 4448, 2991, 2849, 4449, 2992, 2850,
        4038, 4530, 2851, 1324, 2852, 4450, 2993, 2853, 4451, 2994, 2854, 
       4452, 2995, 2855, 4453, 2996, 2856, 4454, 2997, 2857}], 
      Line3DBox[{2859, 4039, 4531, 2858, 4040, 4532, 2860, 4455, 2998, 2861, 
       4456, 2999, 2862, 4457, 3000, 2863, 4458, 3001, 2864, 4459, 3002, 2865,
        4041, 4533, 2866, 4042, 4534, 2867, 1339, 2868, 4460, 3003, 2869, 
       4461, 3004, 2870, 4462, 3005, 2871, 4463, 3006, 2872}], 
      Line3DBox[{2874, 4043, 4535, 2873, 4044, 4536, 2875, 4045, 4537, 2876, 
       4464, 3007, 2877, 4465, 3008, 2878, 4466, 3009, 2879, 4467, 3010, 2880,
        4046, 4538, 2881, 4047, 4539, 2882, 4048, 4540, 2883, 1354, 2884, 
       4468, 3011, 2885, 4469, 3012, 2886, 4470, 3013, 2887}], 
      Line3DBox[{2889, 4049, 4541, 2888, 4050, 4542, 2890, 4051, 4543, 2891, 
       4052, 4544, 2892, 4471, 3014, 2893, 4472, 3015, 2894, 4473, 3016, 2895,
        4053, 4545, 2896, 4054, 4546, 2897, 4055, 4547, 2898, 4056, 4548, 
       2899, 1369, 2900, 4474, 3017, 2901, 4475, 3018, 2902}], 
      Line3DBox[{2904, 4057, 4549, 2903, 4058, 4550, 2905, 4059, 4551, 2906, 
       4060, 4552, 2907, 4061, 4553, 2908, 4476, 3019, 2909, 4477, 3020, 2910,
        4062, 4554, 2911, 4063, 4555, 2912, 4064, 4556, 2913, 4065, 4557, 
       2914, 4066, 4558, 2915, 1384, 2916, 4478, 3021, 2917}], 
      Line3DBox[{2931, 3027, 1403, 2930, 4569, 4077, 2929, 4568, 4076, 2928, 
       4567, 4075, 2927, 4566, 4074, 2926, 4565, 4073, 2925, 4564, 4072, 2924,
        3022, 4479, 2923, 4563, 4071, 2922, 4562, 4070, 2921, 4561, 4069, 
       2920, 4560, 4068, 2919, 4559, 4067, 2918, 1402, 3025, 3026}], 
      Line3DBox[{3029, 3342, 1630, 3028, 4675, 3251, 3030, 4676, 3252, 3031, 
       4677, 3253, 3032, 4678, 3254, 3033, 4679, 3255, 3034, 4680, 3256, 3035,
        4681, 4078, 4570, 3036, 4682, 3257, 3037, 4683, 3258, 3038, 4684, 
       3259, 3039, 4685, 3260, 3040, 4686, 3261, 3041, 4181, 4779, 3262, 
       3343}], Line3DBox[{3043, 4079, 4571, 3042, 1645, 3044, 4687, 3263, 
       3045, 4688, 3264, 3046, 4689, 3265, 3047, 4690, 3266, 3048, 4691, 3267,
        3049, 4692, 4080, 4572, 3050, 4081, 4573, 3051, 4693, 3268, 3052, 
       4694, 3269, 3053, 4695, 3270, 3054, 4696, 3271, 3055, 4697, 3272, 
       3056}], Line3DBox[{3058, 4082, 4574, 3057, 4083, 4575, 3059, 1661, 
       3060, 4698, 3273, 3061, 4699, 3274, 3062, 4700, 3275, 3063, 4701, 3276,
        3064, 4702, 4084, 4576, 3065, 4085, 4577, 3066, 4086, 4578, 3067, 
       4703, 3277, 3068, 4704, 3278, 3069, 4705, 3279, 3070, 4706, 3280, 
       3071}], Line3DBox[{3073, 4087, 4579, 3072, 4088, 4580, 3074, 4089, 
       4581, 3075, 1677, 3076, 4707, 3281, 3077, 4708, 3282, 3078, 4709, 3283,
        3079, 4710, 4090, 4582, 3080, 4091, 4583, 3081, 4092, 4584, 3082, 
       4093, 4585, 3083, 4711, 3284, 3084, 4712, 3285, 3085, 4713, 3286, 
       3086}], Line3DBox[{3088, 4094, 4586, 3087, 4095, 4587, 3089, 4096, 
       4588, 3090, 4097, 4589, 3091, 1693, 3092, 4714, 3287, 3093, 4715, 3288,
        3094, 4716, 4098, 4590, 3095, 4099, 4591, 3096, 4100, 4592, 3097, 
       4101, 4593, 3098, 4102, 4594, 3099, 4717, 3289, 3100, 4718, 3290, 
       3101}], Line3DBox[{3103, 4103, 4595, 3102, 4104, 4596, 3104, 4105, 
       4597, 3105, 4106, 4598, 3106, 4107, 4599, 3107, 1709, 3108, 4719, 3291,
        3109, 4720, 4108, 4600, 3110, 4109, 4601, 3111, 4110, 4602, 3112, 
       4111, 4603, 3113, 4112, 4604, 3114, 4113, 4605, 3115, 4721, 3292, 
       3116}], Line3DBox[{3118, 4114, 4606, 3117, 4115, 4607, 3119, 4116, 
       4608, 3120, 4117, 4609, 3121, 4118, 4610, 3122, 4119, 4611, 3123, 1725,
        3124, 4722, 4120, 4612, 3125, 4121, 4613, 3126, 4122, 4614, 3127, 
       4123, 4615, 3128, 4124, 4616, 3129, 4125, 4617, 3130, 4126, 4618, 
       3131}], Line3DBox[{3135, 4723, 3293, 3133, 4724, 3294, 3137, 4725, 
       3295, 3139, 4726, 3296, 3141, 4727, 3297, 3143, 4728, 3298, 3145, 4729,
        3299, 3147, 4730, 1742, 3149, 4731, 3300, 3151, 4732, 3301, 3153, 
       4733, 3302, 3155, 4734, 3303, 3157, 4735, 3304, 3159, 4736, 3305, 
       3161}], Line3DBox[{3160, 4632, 4139, 3158, 4631, 4138, 3156, 4630, 
       4137, 3154, 4629, 4136, 3152, 4628, 4135, 3150, 4627, 4134, 3148, 4626,
        1741, 3146, 4625, 4133, 3144, 4624, 4132, 3142, 4623, 4131, 3140, 
       4622, 4130, 3138, 4621, 4129, 3136, 4620, 4128, 3132, 4619, 4127, 
       3134}], Line3DBox[{3163, 4140, 4633, 3162, 4737, 3306, 3164, 4738, 
       3307, 3165, 4739, 3308, 3166, 4740, 3309, 3167, 4741, 3310, 3168, 4742,
        3311, 3169, 4743, 4141, 4634, 3170, 1758, 3171, 4744, 3312, 3172, 
       4745, 3313, 3173, 4746, 3314, 3174, 4747, 3315, 3175, 4748, 3316, 
       3176}], Line3DBox[{3178, 4142, 4635, 3177, 4143, 4636, 3179, 4749, 
       3317, 3180, 4750, 3318, 3181, 4751, 3319, 3182, 4752, 3320, 3183, 4753,
        3321, 3184, 4754, 4144, 4637, 3185, 4145, 4638, 3186, 1774, 3187, 
       4755, 3322, 3188, 4756, 3323, 3189, 4757, 3324, 3190, 4758, 3325, 
       3191}], Line3DBox[{3193, 4146, 4639, 3192, 4147, 4640, 3194, 4148, 
       4641, 3195, 4759, 3326, 3196, 4760, 3327, 3197, 4761, 3328, 3198, 4762,
        3329, 3199, 4763, 4149, 4642, 3200, 4150, 4643, 3201, 4151, 4644, 
       3202, 1790, 3203, 4764, 3330, 3204, 4765, 3331, 3205, 4766, 3332, 
       3206}], Line3DBox[{3208, 4152, 4645, 3207, 4153, 4646, 3209, 4154, 
       4647, 3210, 4155, 4648, 3211, 4767, 3333, 3212, 4768, 3334, 3213, 4769,
        3335, 3214, 4770, 4156, 4649, 3215, 4157, 4650, 3216, 4158, 4651, 
       3217, 4159, 4652, 3218, 1806, 3219, 4771, 3336, 3220, 4772, 3337, 
       3221}], Line3DBox[{3223, 4160, 4653, 3222, 4161, 4654, 3224, 4162, 
       4655, 3225, 4163, 4656, 3226, 4164, 4657, 3227, 4773, 3338, 3228, 4774,
        3339, 3229, 4775, 4165, 4658, 3230, 4166, 4659, 3231, 4167, 4660, 
       3232, 4168, 4661, 3233, 4169, 4662, 3234, 1822, 3235, 4776, 3340, 
       3236}], Line3DBox[{3250, 3347, 1842, 3249, 4674, 4180, 3248, 4673, 
       4179, 3247, 4672, 4178, 3246, 4671, 4177, 3245, 4670, 4176, 3244, 4669,
        4175, 4778, 3243, 3341, 4777, 3242, 4668, 4174, 3241, 4667, 4173, 
       3240, 4666, 4172, 3239, 4665, 4171, 3238, 4664, 4170, 3237, 4663, 3345,
        3344, 3346}], 
      Line3DBox[{3349, 3662, 2069, 3348, 4885, 3571, 3350, 4886, 3572, 3351, 
       4887, 3573, 3352, 4888, 3574, 3353, 4889, 3575, 3354, 4890, 3576, 3355,
        4891, 4182, 4780, 3356, 4892, 3577, 3357, 4893, 3578, 3358, 4894, 
       3579, 3359, 4895, 3580, 3360, 4896, 3581, 3361, 4285, 4989, 3582, 
       3663}], Line3DBox[{3363, 4183, 4781, 3362, 2084, 3364, 4897, 3583, 
       3365, 4898, 3584, 3366, 4899, 3585, 3367, 4900, 3586, 3368, 4901, 3587,
        3369, 4902, 4184, 4782, 3370, 4185, 4783, 3371, 4903, 3588, 3372, 
       4904, 3589, 3373, 4905, 3590, 3374, 4906, 3591, 3375, 4907, 3592, 
       3376}], Line3DBox[{3378, 4186, 4784, 3377, 4187, 4785, 3379, 2100, 
       3380, 4908, 3593, 3381, 4909, 3594, 3382, 4910, 3595, 3383, 4911, 3596,
        3384, 4912, 4188, 4786, 3385, 4189, 4787, 3386, 4190, 4788, 3387, 
       4913, 3597, 3388, 4914, 3598, 3389, 4915, 3599, 3390, 4916, 3600, 
       3391}], Line3DBox[{3393, 4191, 4789, 3392, 4192, 4790, 3394, 4193, 
       4791, 3395, 2116, 3396, 4917, 3601, 3397, 4918, 3602, 3398, 4919, 3603,
        3399, 4920, 4194, 4792, 3400, 4195, 4793, 3401, 4196, 4794, 3402, 
       4197, 4795, 3403, 4921, 3604, 3404, 4922, 3605, 3405, 4923, 3606, 
       3406}], Line3DBox[{3408, 4198, 4796, 3407, 4199, 4797, 3409, 4200, 
       4798, 3410, 4201, 4799, 3411, 2132, 3412, 4924, 3607, 3413, 4925, 3608,
        3414, 4926, 4202, 4800, 3415, 4203, 4801, 3416, 4204, 4802, 3417, 
       4205, 4803, 3418, 4206, 4804, 3419, 4927, 3609, 3420, 4928, 3610, 
       3421}], 
      Line3DBox[{3423, 4207, 4805, 3422, 4208, 4806, 3424, 4209, 4807, 3425, 
       4210, 4808, 3426, 4211, 4809, 3427, 2148, 3428, 4929, 3611, 3429, 4930,
        4212, 4810, 3430, 4213, 4811, 3431, 4214, 4812, 3432, 4215, 4813, 
       3433, 4216, 4814, 3434, 4217, 4815, 3435, 4931, 3612, 3436}], 
      Line3DBox[{3438, 4218, 4816, 3437, 4219, 4817, 3439, 4220, 4818, 3440, 
       4221, 4819, 3441, 4222, 4820, 3442, 4223, 4821, 3443, 2164, 3444, 4932,
        4224, 4822, 3445, 4225, 4823, 3446, 4226, 4824, 3447, 4227, 4825, 
       3448, 4228, 4826, 3449, 4229, 4827, 3450, 4230, 4828, 3451}], 
      Line3DBox[{3455, 4933, 3613, 3453, 4934, 3614, 3457, 4935, 3615, 3459, 
       4936, 3616, 3461, 4937, 3617, 3463, 4938, 3618, 3465, 4939, 3619, 3467,
        4940, 2181, 3469, 4941, 3620, 3471, 4942, 3621, 3473, 4943, 3622, 
       3475, 4944, 3623, 3477, 4945, 3624, 3479, 4946, 3625, 3481}], 
      Line3DBox[{3480, 4842, 4243, 3478, 4841, 4242, 3476, 4840, 4241, 3474, 
       4839, 4240, 3472, 4838, 4239, 3470, 4837, 4238, 3468, 4836, 2180, 3466,
        4835, 4237, 3464, 4834, 4236, 3462, 4833, 4235, 3460, 4832, 4234, 
       3458, 4831, 4233, 3456, 4830, 4232, 3452, 4829, 4231, 3454}], 
      Line3DBox[{3483, 4244, 4843, 3482, 4947, 3626, 3484, 4948, 3627, 3485, 
       4949, 3628, 3486, 4950, 3629, 3487, 4951, 3630, 3488, 4952, 3631, 3489,
        4953, 4245, 4844, 3490, 2197, 3491, 4954, 3632, 3492, 4955, 3633, 
       3493, 4956, 3634, 3494, 4957, 3635, 3495, 4958, 3636, 3496}], 
      Line3DBox[{3498, 4246, 4845, 3497, 4247, 4846, 3499, 4959, 3637, 3500, 
       4960, 3638, 3501, 4961, 3639, 3502, 4962, 3640, 3503, 4963, 3641, 3504,
        4964, 4248, 4847, 3505, 4249, 4848, 3506, 2213, 3507, 4965, 3642, 
       3508, 4966, 3643, 3509, 4967, 3644, 3510, 4968, 3645, 3511}], 
      Line3DBox[{3513, 4250, 4849, 3512, 4251, 4850, 3514, 4252, 4851, 3515, 
       4969, 3646, 3516, 4970, 3647, 3517, 4971, 3648, 3518, 4972, 3649, 3519,
        4973, 4253, 4852, 3520, 4254, 4853, 3521, 4255, 4854, 3522, 2229, 
       3523, 4974, 3650, 3524, 4975, 3651, 3525, 4976, 3652, 3526}], 
      Line3DBox[{3528, 4256, 4855, 3527, 4257, 4856, 3529, 4258, 4857, 3530, 
       4259, 4858, 3531, 4977, 3653, 3532, 4978, 3654, 3533, 4979, 3655, 3534,
        4980, 4260, 4859, 3535, 4261, 4860, 3536, 4262, 4861, 3537, 4263, 
       4862, 3538, 2245, 3539, 4981, 3656, 3540, 4982, 3657, 3541}], 
      Line3DBox[{3543, 4264, 4863, 3542, 4265, 4864, 3544, 4266, 4865, 3545, 
       4267, 4866, 3546, 4268, 4867, 3547, 4983, 3658, 3548, 4984, 3659, 3549,
        4985, 4269, 4868, 3550, 4270, 4869, 3551, 4271, 4870, 3552, 4272, 
       4871, 3553, 4273, 4872, 3554, 2261, 3555, 4986, 3660, 3556}], 
      Line3DBox[{3570, 3667, 2281, 3569, 4884, 4284, 3568, 4883, 4283, 3567, 
       4882, 4282, 3566, 4881, 4281, 3565, 4880, 4280, 3564, 4879, 4279, 4988,
        3563, 3661, 4987, 3562, 4878, 4278, 3561, 4877, 4277, 3560, 4876, 
       4276, 3559, 4875, 4275, 3558, 4874, 4274, 3557, 4873, 3665, 3664, 
       3666}], Line3DBox[{3669, 3982, 2508, 3668, 5095, 3891, 3670, 5096, 
       3892, 3671, 5097, 3893, 3672, 5098, 3894, 3673, 5099, 3895, 3674, 5100,
        3896, 3675, 5101, 4286, 4990, 3676, 5102, 3897, 3677, 5103, 3898, 
       3678, 5104, 3899, 3679, 5105, 3900, 3680, 5106, 3901, 3681, 4389, 5199,
        3902, 3983}], 
      Line3DBox[{3683, 4287, 4991, 3682, 2523, 3684, 5107, 3903, 3685, 5108, 
       3904, 3686, 5109, 3905, 3687, 5110, 3906, 3688, 5111, 3907, 3689, 5112,
        4288, 4992, 3690, 4289, 4993, 3691, 5113, 3908, 3692, 5114, 3909, 
       3693, 5115, 3910, 3694, 5116, 3911, 3695, 5117, 3912, 3696}], 
      Line3DBox[{3698, 4290, 4994, 3697, 4291, 4995, 3699, 2539, 3700, 5118, 
       3913, 3701, 5119, 3914, 3702, 5120, 3915, 3703, 5121, 3916, 3704, 5122,
        4292, 4996, 3705, 4293, 4997, 3706, 4294, 4998, 3707, 5123, 3917, 
       3708, 5124, 3918, 3709, 5125, 3919, 3710, 5126, 3920, 3711}], 
      Line3DBox[{3713, 4295, 4999, 3712, 4296, 5000, 3714, 4297, 5001, 3715, 
       2555, 3716, 5127, 3921, 3717, 5128, 3922, 3718, 5129, 3923, 3719, 5130,
        4298, 5002, 3720, 4299, 5003, 3721, 4300, 5004, 3722, 4301, 5005, 
       3723, 5131, 3924, 3724, 5132, 3925, 3725, 5133, 3926, 3726}], 
      Line3DBox[{3728, 4302, 5006, 3727, 4303, 5007, 3729, 4304, 5008, 3730, 
       4305, 5009, 3731, 2571, 3732, 5134, 3927, 3733, 5135, 3928, 3734, 5136,
        4306, 5010, 3735, 4307, 5011, 3736, 4308, 5012, 3737, 4309, 5013, 
       3738, 4310, 5014, 3739, 5137, 3929, 3740, 5138, 3930, 3741}], 
      Line3DBox[{3743, 4311, 5015, 3742, 4312, 5016, 3744, 4313, 5017, 3745, 
       4314, 5018, 3746, 4315, 5019, 3747, 2587, 3748, 5139, 3931, 3749, 5140,
        4316, 5020, 3750, 4317, 5021, 3751, 4318, 5022, 3752, 4319, 5023, 
       3753, 4320, 5024, 3754, 4321, 5025, 3755, 5141, 3932, 3756}], 
      Line3DBox[{3758, 4322, 5026, 3757, 4323, 5027, 3759, 4324, 5028, 3760, 
       4325, 5029, 3761, 4326, 5030, 3762, 4327, 5031, 3763, 2603, 3764, 5142,
        4328, 5032, 3765, 4329, 5033, 3766, 4330, 5034, 3767, 4331, 5035, 
       3768, 4332, 5036, 3769, 4333, 5037, 3770, 4334, 5038, 3771}], 
      Line3DBox[{3775, 5143, 3933, 3773, 5144, 3934, 3777, 5145, 3935, 3779, 
       5146, 3936, 3781, 5147, 3937, 3783, 5148, 3938, 3785, 5149, 3939, 3787,
        5150, 2620, 3789, 5151, 3940, 3791, 5152, 3941, 3793, 5153, 3942, 
       3795, 5154, 3943, 3797, 5155, 3944, 3799, 5156, 3945, 3801}], 
      Line3DBox[{3800, 5052, 4347, 3798, 5051, 4346, 3796, 5050, 4345, 3794, 
       5049, 4344, 3792, 5048, 4343, 3790, 5047, 4342, 3788, 5046, 2619, 3786,
        5045, 4341, 3784, 5044, 4340, 3782, 5043, 4339, 3780, 5042, 4338, 
       3778, 5041, 4337, 3776, 5040, 4336, 3772, 5039, 4335, 3774}], 
      Line3DBox[{3803, 4348, 5053, 3802, 5157, 3946, 3804, 5158, 3947, 3805, 
       5159, 3948, 3806, 5160, 3949, 3807, 5161, 3950, 3808, 5162, 3951, 3809,
        5163, 4349, 5054, 3810, 2636, 3811, 5164, 3952, 3812, 5165, 3953, 
       3813, 5166, 3954, 3814, 5167, 3955, 3815, 5168, 3956, 3816}], 
      Line3DBox[{3818, 4350, 5055, 3817, 4351, 5056, 3819, 5169, 3957, 3820, 
       5170, 3958, 3821, 5171, 3959, 3822, 5172, 3960, 3823, 5173, 3961, 3824,
        5174, 4352, 5057, 3825, 4353, 5058, 3826, 2652, 3827, 5175, 3962, 
       3828, 5176, 3963, 3829, 5177, 3964, 3830, 5178, 3965, 3831}], 
      Line3DBox[{3833, 4354, 5059, 3832, 4355, 5060, 3834, 4356, 5061, 3835, 
       5179, 3966, 3836, 5180, 3967, 3837, 5181, 3968, 3838, 5182, 3969, 3839,
        5183, 4357, 5062, 3840, 4358, 5063, 3841, 4359, 5064, 3842, 2668, 
       3843, 5184, 3970, 3844, 5185, 3971, 3845, 5186, 3972, 3846}], 
      Line3DBox[{3848, 4360, 5065, 3847, 4361, 5066, 3849, 4362, 5067, 3850, 
       4363, 5068, 3851, 5187, 3973, 3852, 5188, 3974, 3853, 5189, 3975, 3854,
        5190, 4364, 5069, 3855, 4365, 5070, 3856, 4366, 5071, 3857, 4367, 
       5072, 3858, 2684, 3859, 5191, 3976, 3860, 5192, 3977, 3861}], 
      Line3DBox[{3863, 4368, 5073, 3862, 4369, 5074, 3864, 4370, 5075, 3865, 
       4371, 5076, 3866, 4372, 5077, 3867, 5193, 3978, 3868, 5194, 3979, 3869,
        5195, 4373, 5078, 3870, 4374, 5079, 3871, 4375, 5080, 3872, 4376, 
       5081, 3873, 4377, 5082, 3874, 2700, 3875, 5196, 3980, 3876}], 
      Line3DBox[{3890, 3987, 2720, 3889, 5094, 4388, 3888, 5093, 4387, 3887, 
       5092, 4386, 3886, 5091, 4385, 3885, 5090, 4384, 3884, 5089, 4383, 5198,
        3883, 3981, 5197, 3882, 5088, 4382, 3881, 5087, 4381, 3880, 5086, 
       4380, 3879, 5085, 4379, 3878, 5084, 4378, 3877, 5083, 3985, 3984, 
       3986}]}, {
      Line3DBox[{998, 1206, 4390, 999, 1219, 1025, 4485, 1233, 1039, 4490, 
       1247, 1053, 4497, 1261, 1067, 4506, 1275, 1081, 4517, 1289, 1095, 1303,
        4432, 1109, 1317, 4444, 1123, 4532, 1331, 1137, 4536, 1345, 1151, 
       4542, 1359, 1165, 4550, 1373, 1179, 4559, 1387, 1193}], 
      Line3DBox[{1000, 1207, 4391, 1001, 1220, 4401, 1026, 1234, 1040, 4491, 
       1248, 1054, 4498, 1262, 1068, 4507, 1276, 1082, 4518, 1290, 1096, 1304,
        4433, 1110, 1318, 4445, 1124, 1332, 4455, 1138, 4537, 1346, 1152, 
       4543, 1360, 1166, 4551, 1374, 1180, 4560, 1388, 1194}], 
      Line3DBox[{1002, 1208, 4392, 1003, 1221, 4402, 1027, 1235, 4411, 1041, 
       1249, 1055, 4499, 1263, 1069, 4508, 1277, 1083, 4519, 1291, 1097, 1305,
        4434, 1111, 1319, 4446, 1125, 1333, 4456, 1139, 1347, 4464, 1153, 
       4544, 1361, 1167, 4552, 1375, 1181, 4561, 1389, 1195}], 
      Line3DBox[{1004, 1209, 4393, 1005, 1222, 4403, 1028, 1236, 4412, 1042, 
       1250, 4419, 1056, 1264, 1070, 4509, 1278, 1084, 4520, 1292, 1098, 1306,
        4435, 1112, 1320, 4447, 1126, 1334, 4457, 1140, 1348, 4465, 1154, 
       1362, 4471, 1168, 4553, 1376, 1182, 4562, 1390, 1196}], 
      Line3DBox[{1006, 1210, 4394, 1007, 1223, 4404, 1029, 1237, 4413, 1043, 
       1251, 4420, 1057, 1265, 4425, 1071, 1279, 1085, 4521, 1293, 1099, 1307,
        4436, 1113, 1321, 4448, 1127, 1335, 4458, 1141, 1349, 4466, 1155, 
       1363, 4472, 1169, 1377, 4476, 1183, 4563, 1391, 1197}], 
      Line3DBox[{1008, 1211, 4395, 1009, 1224, 4405, 1030, 1238, 4414, 1044, 
       1252, 4421, 1058, 1266, 4426, 1072, 1280, 4429, 1086, 1294, 1100, 1308,
        4437, 1114, 1322, 4449, 1128, 1336, 4459, 1142, 1350, 4467, 1156, 
       1364, 4473, 1170, 1378, 4477, 1184, 1392, 4479, 1198}], 
      Line3DBox[{8, 2731, 23, 2745, 38, 2760, 53, 2775, 68, 2790, 83, 2805, 
       98, 2820, 113, 2835, 128, 2850, 143, 2865, 158, 2880, 173, 2895, 188, 
       2910, 203, 2924, 218}], 
      Line3DBox[{1010, 4480, 1212, 1011, 4482, 1225, 1031, 4486, 1239, 1045, 
       4492, 1253, 1059, 4500, 1267, 1073, 4510, 1281, 1087, 4522, 1295, 1101,
        1309, 1115, 4530, 1323, 1129, 4533, 1337, 1143, 4538, 1351, 1157, 
       4545, 1365, 1171, 4554, 1379, 1185, 4564, 1393, 1199}], 
      Line3DBox[{1012, 1213, 4396, 1013, 4483, 1226, 1032, 4487, 1240, 1046, 
       4493, 1254, 1060, 4501, 1268, 1074, 4511, 1282, 1088, 4523, 1296, 1102,
        1310, 4438, 1116, 1324, 1130, 4534, 1338, 1144, 4539, 1352, 1158, 
       4546, 1366, 1172, 4555, 1380, 1186, 4565, 1394, 1200}], 
      Line3DBox[{1014, 1214, 4397, 1015, 1227, 4406, 1033, 4488, 1241, 1047, 
       4494, 1255, 1061, 4502, 1269, 1075, 4512, 1283, 1089, 4524, 1297, 1103,
        1311, 4439, 1117, 1325, 4450, 1131, 1339, 1145, 4540, 1353, 1159, 
       4547, 1367, 1173, 4556, 1381, 1187, 4566, 1395, 1201}], 
      Line3DBox[{1016, 1215, 4398, 1017, 1228, 4407, 1034, 1242, 4415, 1048, 
       4495, 1256, 1062, 4503, 1270, 1076, 4513, 1284, 1090, 4525, 1298, 1104,
        1312, 4440, 1118, 1326, 4451, 1132, 1340, 4460, 1146, 1354, 1160, 
       4548, 1368, 1174, 4557, 1382, 1188, 4567, 1396, 1202}], 
      Line3DBox[{1018, 1216, 4399, 1019, 1229, 4408, 1035, 1243, 4416, 1049, 
       1257, 4422, 1063, 4504, 1271, 1077, 4514, 1285, 1091, 4526, 1299, 1105,
        1313, 4441, 1119, 1327, 4452, 1133, 1341, 4461, 1147, 1355, 4468, 
       1161, 1369, 1175, 4558, 1383, 1189, 4568, 1397, 1203}], 
      Line3DBox[{1020, 1217, 4400, 1021, 1230, 4409, 1036, 1244, 4417, 1050, 
       1258, 4423, 1064, 1272, 4427, 1078, 4515, 1286, 1092, 4527, 1300, 1106,
        1314, 4442, 1120, 1328, 4453, 1134, 1342, 4462, 1148, 1356, 4469, 
       1162, 1370, 4474, 1176, 1384, 1190, 4569, 1398, 1204}], 
      Line3DBox[{1022, 1400, 1401, 1023, 1231, 4410, 1037, 1245, 4418, 1051, 
       1259, 4424, 1065, 1273, 4428, 1079, 1287, 4430, 1093, 4528, 1301, 1107,
        1315, 4443, 1121, 1329, 4454, 1135, 1343, 4463, 1149, 1357, 4470, 
       1163, 1371, 4475, 1177, 1385, 4478, 1191, 1403, 1404, 1405}], 
      Line3DBox[{1192, 1386, 1402, 1178, 1372, 4549, 1164, 1358, 4541, 1150, 
       1344, 4535, 1136, 1330, 4531, 1122, 1316, 4529, 1108, 4431, 1302, 1094,
        1288, 4516, 1080, 1274, 4505, 1066, 1260, 4496, 1052, 1246, 4489, 
       1038, 1232, 4484, 1024, 1218, 4481, 997, 1205, 1399, 1406}], 
      Line3DBox[{1408, 1631, 4675, 1409, 1645, 1437, 4575, 1660, 1452, 4580, 
       1675, 1467, 4587, 1690, 1482, 4596, 1705, 1497, 4607, 1720, 1512, 4620,
        1735, 4724, 1527, 1750, 4737, 1542, 4636, 1765, 1557, 4640, 1780, 
       1572, 4646, 1795, 1587, 4654, 1810, 1602, 4664, 1825, 1617}], 
      Line3DBox[{1410, 1632, 4676, 1411, 1646, 4687, 1438, 1661, 1453, 4581, 
       1676, 1468, 4588, 1691, 1483, 4597, 1706, 1498, 4608, 1721, 1513, 4621,
        1736, 4725, 1528, 1751, 4738, 1543, 1766, 4749, 1558, 4641, 1781, 
       1573, 4647, 1796, 1588, 4655, 1811, 1603, 4665, 1826, 1618}], 
      Line3DBox[{1412, 1633, 4677, 1413, 1647, 4688, 1439, 1662, 4698, 1454, 
       1677, 1469, 4589, 1692, 1484, 4598, 1707, 1499, 4609, 1722, 1514, 4622,
        1737, 4726, 1529, 1752, 4739, 1544, 1767, 4750, 1559, 1782, 4759, 
       1574, 4648, 1797, 1589, 4656, 1812, 1604, 4666, 1827, 1619}], 
      Line3DBox[{1414, 1634, 4678, 1415, 1648, 4689, 1440, 1663, 4699, 1455, 
       1678, 4707, 1470, 1693, 1485, 4599, 1708, 1500, 4610, 1723, 1515, 4623,
        1738, 4727, 1530, 1753, 4740, 1545, 1768, 4751, 1560, 1783, 4760, 
       1575, 1798, 4767, 1590, 4657, 1813, 1605, 4667, 1828, 1620}], 
      Line3DBox[{1416, 1635, 4679, 1417, 1649, 4690, 1441, 1664, 4700, 1456, 
       1679, 4708, 1471, 1694, 4714, 1486, 1709, 1501, 4611, 1724, 1516, 4624,
        1739, 4728, 1531, 1754, 4741, 1546, 1769, 4752, 1561, 1784, 4761, 
       1576, 1799, 4768, 1591, 1814, 4773, 1606, 4668, 1829, 1621}], 
      Line3DBox[{1418, 1636, 4680, 1419, 1650, 4691, 1442, 1665, 4701, 1457, 
       1680, 4709, 1472, 1695, 4715, 1487, 1710, 4719, 1502, 1725, 1517, 4625,
        1740, 4729, 1532, 1755, 4742, 1547, 1770, 4753, 1562, 1785, 4762, 
       1577, 1800, 4769, 1592, 1815, 4774, 1607, 1830, 4777, 1622}], 
      Line3DBox[{1420, 1637, 4681, 1422, 1651, 4692, 1443, 1666, 4702, 1458, 
       1681, 4710, 1473, 1696, 4716, 1488, 1711, 4720, 1503, 1726, 4722, 1518,
        1741, 4730, 1533, 1756, 4743, 1548, 1771, 4754, 1563, 1786, 4763, 
       1578, 1801, 4770, 1593, 1816, 4775, 1608, 1831, 4778, 1623}], 
      Line3DBox[{1424, 1639, 4682, 1425, 4573, 1653, 1445, 4577, 1668, 1460, 
       4583, 1683, 1475, 4591, 1698, 1490, 4601, 1713, 1505, 4613, 1728, 1520,
        4627, 1743, 4731, 1535, 1758, 1550, 4638, 1773, 1565, 4643, 1788, 
       1580, 4650, 1803, 1595, 4659, 1818, 1610, 4670, 1833, 1625}], 
      Line3DBox[{1426, 1640, 4683, 1427, 1654, 4693, 1446, 4578, 1669, 1461, 
       4584, 1684, 1476, 4592, 1699, 1491, 4602, 1714, 1506, 4614, 1729, 1521,
        4628, 1744, 4732, 1536, 1759, 4744, 1551, 1774, 1566, 4644, 1789, 
       1581, 4651, 1804, 1596, 4660, 1819, 1611, 4671, 1834, 1626}], 
      Line3DBox[{1428, 1641, 4684, 1429, 1655, 4694, 1447, 1670, 4703, 1462, 
       4585, 1685, 1477, 4593, 1700, 1492, 4603, 1715, 1507, 4615, 1730, 1522,
        4629, 1745, 4733, 1537, 1760, 4745, 1552, 1775, 4755, 1567, 1790, 
       1582, 4652, 1805, 1597, 4661, 1820, 1612, 4672, 1835, 1627}], 
      Line3DBox[{1430, 1642, 4685, 1431, 1656, 4695, 1448, 1671, 4704, 1463, 
       1686, 4711, 1478, 4594, 1701, 1493, 4604, 1716, 1508, 4616, 1731, 1523,
        4630, 1746, 4734, 1538, 1761, 4746, 1553, 1776, 4756, 1568, 1791, 
       4764, 1583, 1806, 1598, 4662, 1821, 1613, 4673, 1836, 1628}], 
      Line3DBox[{1432, 1643, 4686, 1433, 1657, 4696, 1449, 1672, 4705, 1464, 
       1687, 4712, 1479, 1702, 4717, 1494, 4605, 1717, 1509, 4617, 1732, 1524,
        4631, 1747, 4735, 1539, 1762, 4747, 1554, 1777, 4757, 1569, 1792, 
       4765, 1584, 1807, 4771, 1599, 1822, 1614, 4674, 1837, 1629}], 
      Line3DBox[{1434, 1839, 1840, 4779, 1435, 1658, 4697, 1450, 1673, 4706, 
       1465, 1688, 4713, 1480, 1703, 4718, 1495, 1718, 4721, 1510, 4618, 1733,
        1525, 4632, 1748, 4736, 1540, 1763, 4748, 1555, 1778, 4758, 1570, 
       1793, 4766, 1585, 1808, 4772, 1600, 1823, 4776, 1615, 1842, 1843, 
       1844}], Line3DBox[{1616, 1824, 4663, 1841, 1601, 1809, 4653, 1586, 
       1794, 4645, 1571, 1779, 4639, 1556, 1764, 4635, 1541, 1749, 4633, 1526,
        4723, 1734, 4619, 1511, 1719, 4606, 1496, 1704, 4595, 1481, 1689, 
       4586, 1466, 1674, 4579, 1451, 1659, 4574, 1436, 1644, 4571, 1407, 1630,
        1838, 1845}], 
      Line3DBox[{1624, 1832, 4669, 1609, 1817, 4658, 1594, 1802, 4649, 1579, 
       1787, 4642, 1564, 1772, 4637, 1549, 1757, 4634, 1534, 1742, 4626, 1519,
        1727, 4612, 1504, 1712, 4600, 1489, 1697, 4590, 1474, 1682, 4582, 
       1459, 1667, 4576, 1444, 1652, 4572, 1423, 1638, 4570, 1421}], 
      Line3DBox[{1847, 2070, 4885, 1848, 2084, 1876, 4785, 2099, 1891, 4790, 
       2114, 1906, 4797, 2129, 1921, 4806, 2144, 1936, 4817, 2159, 1951, 4830,
        2174, 4934, 1966, 2189, 4947, 1981, 4846, 2204, 1996, 4850, 2219, 
       2011, 4856, 2234, 2026, 4864, 2249, 2041, 4874, 2264, 2056}], 
      Line3DBox[{1849, 2071, 4886, 1850, 2085, 4897, 1877, 2100, 1892, 4791, 
       2115, 1907, 4798, 2130, 1922, 4807, 2145, 1937, 4818, 2160, 1952, 4831,
        2175, 4935, 1967, 2190, 4948, 1982, 2205, 4959, 1997, 4851, 2220, 
       2012, 4857, 2235, 2027, 4865, 2250, 2042, 4875, 2265, 2057}], 
      Line3DBox[{1851, 2072, 4887, 1852, 2086, 4898, 1878, 2101, 4908, 1893, 
       2116, 1908, 4799, 2131, 1923, 4808, 2146, 1938, 4819, 2161, 1953, 4832,
        2176, 4936, 1968, 2191, 4949, 1983, 2206, 4960, 1998, 2221, 4969, 
       2013, 4858, 2236, 2028, 4866, 2251, 2043, 4876, 2266, 2058}], 
      Line3DBox[{1853, 2073, 4888, 1854, 2087, 4899, 1879, 2102, 4909, 1894, 
       2117, 4917, 1909, 2132, 1924, 4809, 2147, 1939, 4820, 2162, 1954, 4833,
        2177, 4937, 1969, 2192, 4950, 1984, 2207, 4961, 1999, 2222, 4970, 
       2014, 2237, 4977, 2029, 4867, 2252, 2044, 4877, 2267, 2059}], 
      Line3DBox[{1855, 2074, 4889, 1856, 2088, 4900, 1880, 2103, 4910, 1895, 
       2118, 4918, 1910, 2133, 4924, 1925, 2148, 1940, 4821, 2163, 1955, 4834,
        2178, 4938, 1970, 2193, 4951, 1985, 2208, 4962, 2000, 2223, 4971, 
       2015, 2238, 4978, 2030, 2253, 4983, 2045, 4878, 2268, 2060}], 
      Line3DBox[{1857, 2075, 4890, 1858, 2089, 4901, 1881, 2104, 4911, 1896, 
       2119, 4919, 1911, 2134, 4925, 1926, 2149, 4929, 1941, 2164, 1956, 4835,
        2179, 4939, 1971, 2194, 4952, 1986, 2209, 4963, 2001, 2224, 4972, 
       2016, 2239, 4979, 2031, 2254, 4984, 2046, 2269, 4987, 2061}], 
      Line3DBox[{1859, 2076, 4891, 1861, 2090, 4902, 1882, 2105, 4912, 1897, 
       2120, 4920, 1912, 2135, 4926, 1927, 2150, 4930, 1942, 2165, 4932, 1957,
        2180, 4940, 1972, 2195, 4953, 1987, 2210, 4964, 2002, 2225, 4973, 
       2017, 2240, 4980, 2032, 2255, 4985, 2047, 2270, 4988, 2062}], 
      Line3DBox[{1863, 2078, 4892, 1864, 4783, 2092, 1884, 4787, 2107, 1899, 
       4793, 2122, 1914, 4801, 2137, 1929, 4811, 2152, 1944, 4823, 2167, 1959,
        4837, 2182, 4941, 1974, 2197, 1989, 4848, 2212, 2004, 4853, 2227, 
       2019, 4860, 2242, 2034, 4869, 2257, 2049, 4880, 2272, 2064}], 
      Line3DBox[{1865, 2079, 4893, 1866, 2093, 4903, 1885, 4788, 2108, 1900, 
       4794, 2123, 1915, 4802, 2138, 1930, 4812, 2153, 1945, 4824, 2168, 1960,
        4838, 2183, 4942, 1975, 2198, 4954, 1990, 2213, 2005, 4854, 2228, 
       2020, 4861, 2243, 2035, 4870, 2258, 2050, 4881, 2273, 2065}], 
      Line3DBox[{1867, 2080, 4894, 1868, 2094, 4904, 1886, 2109, 4913, 1901, 
       4795, 2124, 1916, 4803, 2139, 1931, 4813, 2154, 1946, 4825, 2169, 1961,
        4839, 2184, 4943, 1976, 2199, 4955, 1991, 2214, 4965, 2006, 2229, 
       2021, 4862, 2244, 2036, 4871, 2259, 2051, 4882, 2274, 2066}], 
      Line3DBox[{1869, 2081, 4895, 1870, 2095, 4905, 1887, 2110, 4914, 1902, 
       2125, 4921, 1917, 4804, 2140, 1932, 4814, 2155, 1947, 4826, 2170, 1962,
        4840, 2185, 4944, 1977, 2200, 4956, 1992, 2215, 4966, 2007, 2230, 
       4974, 2022, 2245, 2037, 4872, 2260, 2052, 4883, 2275, 2067}], 
      Line3DBox[{1871, 2082, 4896, 1872, 2096, 4906, 1888, 2111, 4915, 1903, 
       2126, 4922, 1918, 2141, 4927, 1933, 4815, 2156, 1948, 4827, 2171, 1963,
        4841, 2186, 4945, 1978, 2201, 4957, 1993, 2216, 4967, 2008, 2231, 
       4975, 2023, 2246, 4981, 2038, 2261, 2053, 4884, 2276, 2068}], 
      Line3DBox[{1873, 2278, 2279, 4989, 1874, 2097, 4907, 1889, 2112, 4916, 
       1904, 2127, 4923, 1919, 2142, 4928, 1934, 2157, 4931, 1949, 4828, 2172,
        1964, 4842, 2187, 4946, 1979, 2202, 4958, 1994, 2217, 4968, 2009, 
       2232, 4976, 2024, 2247, 4982, 2039, 2262, 4986, 2054, 2281, 2282, 
       2283}], Line3DBox[{2055, 2263, 4873, 2280, 2040, 2248, 4863, 2025, 
       2233, 4855, 2010, 2218, 4849, 1995, 2203, 4845, 1980, 2188, 4843, 1965,
        4933, 2173, 4829, 1950, 2158, 4816, 1935, 2143, 4805, 1920, 2128, 
       4796, 1905, 2113, 4789, 1890, 2098, 4784, 1875, 2083, 4781, 1846, 2069,
        2277, 2284}], 
      Line3DBox[{2063, 2271, 4879, 2048, 2256, 4868, 2033, 2241, 4859, 2018, 
       2226, 4852, 2003, 2211, 4847, 1988, 2196, 4844, 1973, 2181, 4836, 1958,
        2166, 4822, 1943, 2151, 4810, 1928, 2136, 4800, 1913, 2121, 4792, 
       1898, 2106, 4786, 1883, 2091, 4782, 1862, 2077, 4780, 1860}], 
      Line3DBox[{2286, 2509, 5095, 2287, 2523, 2315, 4995, 2538, 2330, 5000, 
       2553, 2345, 5007, 2568, 2360, 5016, 2583, 2375, 5027, 2598, 2390, 5040,
        2613, 5144, 2405, 2628, 5157, 2420, 5056, 2643, 2435, 5060, 2658, 
       2450, 5066, 2673, 2465, 5074, 2688, 2480, 5084, 2703, 2495}], 
      Line3DBox[{2288, 2510, 5096, 2289, 2524, 5107, 2316, 2539, 2331, 5001, 
       2554, 2346, 5008, 2569, 2361, 5017, 2584, 2376, 5028, 2599, 2391, 5041,
        2614, 5145, 2406, 2629, 5158, 2421, 2644, 5169, 2436, 5061, 2659, 
       2451, 5067, 2674, 2466, 5075, 2689, 2481, 5085, 2704, 2496}], 
      Line3DBox[{2290, 2511, 5097, 2291, 2525, 5108, 2317, 2540, 5118, 2332, 
       2555, 2347, 5009, 2570, 2362, 5018, 2585, 2377, 5029, 2600, 2392, 5042,
        2615, 5146, 2407, 2630, 5159, 2422, 2645, 5170, 2437, 2660, 5179, 
       2452, 5068, 2675, 2467, 5076, 2690, 2482, 5086, 2705, 2497}], 
      Line3DBox[{2292, 2512, 5098, 2293, 2526, 5109, 2318, 2541, 5119, 2333, 
       2556, 5127, 2348, 2571, 2363, 5019, 2586, 2378, 5030, 2601, 2393, 5043,
        2616, 5147, 2408, 2631, 5160, 2423, 2646, 5171, 2438, 2661, 5180, 
       2453, 2676, 5187, 2468, 5077, 2691, 2483, 5087, 2706, 2498}], 
      Line3DBox[{2294, 2513, 5099, 2295, 2527, 5110, 2319, 2542, 5120, 2334, 
       2557, 5128, 2349, 2572, 5134, 2364, 2587, 2379, 5031, 2602, 2394, 5044,
        2617, 5148, 2409, 2632, 5161, 2424, 2647, 5172, 2439, 2662, 5181, 
       2454, 2677, 5188, 2469, 2692, 5193, 2484, 5088, 2707, 2499}], 
      Line3DBox[{2296, 2514, 5100, 2297, 2528, 5111, 2320, 2543, 5121, 2335, 
       2558, 5129, 2350, 2573, 5135, 2365, 2588, 5139, 2380, 2603, 2395, 5045,
        2618, 5149, 2410, 2633, 5162, 2425, 2648, 5173, 2440, 2663, 5182, 
       2455, 2678, 5189, 2470, 2693, 5194, 2485, 2708, 5197, 2500}], 
      Line3DBox[{2298, 2515, 5101, 2300, 2529, 5112, 2321, 2544, 5122, 2336, 
       2559, 5130, 2351, 2574, 5136, 2366, 2589, 5140, 2381, 2604, 5142, 2396,
        2619, 5150, 2411, 2634, 5163, 2426, 2649, 5174, 2441, 2664, 5183, 
       2456, 2679, 5190, 2471, 2694, 5195, 2486, 2709, 5198, 2501}], 
      Line3DBox[{2302, 2517, 5102, 2303, 4993, 2531, 2323, 4997, 2546, 2338, 
       5003, 2561, 2353, 5011, 2576, 2368, 5021, 2591, 2383, 5033, 2606, 2398,
        5047, 2621, 5151, 2413, 2636, 2428, 5058, 2651, 2443, 5063, 2666, 
       2458, 5070, 2681, 2473, 5079, 2696, 2488, 5090, 2711, 2503}], 
      Line3DBox[{2304, 2518, 5103, 2305, 2532, 5113, 2324, 4998, 2547, 2339, 
       5004, 2562, 2354, 5012, 2577, 2369, 5022, 2592, 2384, 5034, 2607, 2399,
        5048, 2622, 5152, 2414, 2637, 5164, 2429, 2652, 2444, 5064, 2667, 
       2459, 5071, 2682, 2474, 5080, 2697, 2489, 5091, 2712, 2504}], 
      Line3DBox[{2306, 2519, 5104, 2307, 2533, 5114, 2325, 2548, 5123, 2340, 
       5005, 2563, 2355, 5013, 2578, 2370, 5023, 2593, 2385, 5035, 2608, 2400,
        5049, 2623, 5153, 2415, 2638, 5165, 2430, 2653, 5175, 2445, 2668, 
       2460, 5072, 2683, 2475, 5081, 2698, 2490, 5092, 2713, 2505}], 
      Line3DBox[{2308, 2520, 5105, 2309, 2534, 5115, 2326, 2549, 5124, 2341, 
       2564, 5131, 2356, 5014, 2579, 2371, 5024, 2594, 2386, 5036, 2609, 2401,
        5050, 2624, 5154, 2416, 2639, 5166, 2431, 2654, 5176, 2446, 2669, 
       5184, 2461, 2684, 2476, 5082, 2699, 2491, 5093, 2714, 2506}], 
      Line3DBox[{2310, 2521, 5106, 2311, 2535, 5116, 2327, 2550, 5125, 2342, 
       2565, 5132, 2357, 2580, 5137, 2372, 5025, 2595, 2387, 5037, 2610, 2402,
        5051, 2625, 5155, 2417, 2640, 5167, 2432, 2655, 5177, 2447, 2670, 
       5185, 2462, 2685, 5191, 2477, 2700, 2492, 5094, 2715, 2507}], 
      Line3DBox[{2312, 2717, 2718, 5199, 2313, 2536, 5117, 2328, 2551, 5126, 
       2343, 2566, 5133, 2358, 2581, 5138, 2373, 2596, 5141, 2388, 5038, 2611,
        2403, 5052, 2626, 5156, 2418, 2641, 5168, 2433, 2656, 5178, 2448, 
       2671, 5186, 2463, 2686, 5192, 2478, 2701, 5196, 2493, 2720, 2721, 
       2722}], Line3DBox[{2494, 2702, 5083, 2719, 2479, 2687, 5073, 2464, 
       2672, 5065, 2449, 2657, 5059, 2434, 2642, 5055, 2419, 2627, 5053, 2404,
        5143, 2612, 5039, 2389, 2597, 5026, 2374, 2582, 5015, 2359, 2567, 
       5006, 2344, 2552, 4999, 2329, 2537, 4994, 2314, 2522, 4991, 2285, 2508,
        2716, 2723}], 
      Line3DBox[{2502, 2710, 5089, 2487, 2695, 5078, 2472, 2680, 5069, 2457, 
       2665, 5062, 2442, 2650, 5057, 2427, 2635, 5054, 2412, 2620, 5046, 2397,
        2605, 5032, 2382, 2590, 5020, 2367, 2575, 5010, 2352, 2560, 5002, 
       2337, 2545, 4996, 2322, 2530, 4992, 2301, 2516, 4990, 
       2299}]}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}, {}}},
   VertexNormals->CompressedData["
1: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
w+6flOiReL7IIbnGzWqlRHGd4l89WPnvKoQoYeFCvdYSEkmvINoy5odEh59L
x508o98DOEjjv4ze92QofQ6Y/X+ezjvkMJ2HmPxzlM77JImuA0z+t6XrLjGm
6zBTfy1p3EOkOIjBP0NoPsPEeZCU19iOmjKe4pPE/FN06x2KX+78ey2GC4Ph
G8XniYPOHt4Jit+f/3stgsZBDeGljUWk0aLbCIm+8vLvtQhu/30/RdL3lcG8
r8gPzkeyrhWRg88DciT64r2/1yIY/vd8FBGDyS66En0Xz4sNraeT7RxdvScd
d5JKvwcm/vPofU+60OeA2f/d6bxDutB5iMk/VnTeJ5vpOsDk/3S67pKNdB1m
6m8ljXuINo2DGPwzi8adRIpDGfxZ5FRhRfFDJs5PJxf+5Ynen48EUfycmC9d
k+5A8fWIbM9eEp5OzK4ojMvIJ6XNCvISvSReZY2SRCfJO/je6nBkPnn6d02C
9xff5kl0qgH6fRwrPPLJkzFKSyV64cJFKqclOuGRoL5rT+/KJ4ldk+0leu2M
Wbt1JTpts2b/lQLIJ2/pFRp5A7pJdPK5dNzJYPo9MPG/Se97so0+B8z+16Dz
DrlO5yEm/6jQeZ9cpusAk/+V6bpLztF1mKm/y2jcI8U/1Qz+6UHjTtJG41AG
f26hcT+R8gAG/+855J9B8W0mzpNezJsn4d2TVq2cln4kh4jnzLmy43cw3HGw
PyLRPcT3noQMccwh6zR+LJfoT2opOsYS3cl061OZc1tzyI45CYWza5PAZdj9
mxLdrzw0THeCUQ6BxwlbJPrr52NuzRLdNSmqNm7l2BzymF6Z85JKx50cpN8D
E/8j9L4n4fQ5YPa/EZ13iD+dh5j8Q+i8T7LoOsDk/yF03SXSOgyXpPWXR+Me
IkvjIAb/NNG4k/SicSiDPzNp3E+8aB7A4H9zmncRKQ9j6qmiVccqNfVMJs6r
hMc8JDrGuj96HlNVM8nqv2sI9NW/qyfRkWalgPelHplE9+8aD9vPx2RIdLwL
ucadQCaTGGt9dJPoqfi+fvou72MvyCDt9AqPLnQfJtGx59NxJyX0e2Din0Lv
e9KfPgfM/jen8w4JovMQk394dN4n5XQdYPK/El13yQW6DjP1ty+Ne8g6Ggcx
+GctjTuJNo1DGfx5nsb95DnNAxj8f4HmXUTKw8Bems/dad5LpDyY4b8Lxn+Z
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
PtEpdt/TFfG8KrvvyeBPMbvvyew3rOtY57G+Yx3CuoT1KIvd92Ty53V235M5
713YfU9mf2IcMa4OnL4nrj7Svqcuu+9JsO+5j933JNj39GX3PQn2PWXYfU+C
fc9Qdt+TYN/zGLuP74r7R47dx2fwfDO7j0+4fXzETYiX1rH7+ITbx8e6hPXo
MruPz+RPPOd47vG8D2X38Zn9qc/u4xNuHx9XxLdm7D4+U78GsPv4BPv459l9
fMLt46OfC31cl9i+FFfkm3VsXwrh+lIQ12N83rN9Kcx+02b7UgjXl4J1Huv7
ObYvhalHm9m+FCZ/4jnHc4/n/TLbl8Lsz21sXwoTz75sXwpBX8pbti+FoC+l
E9uXQtCX4sj2pRD0pSxm+6xckb8jz0TeiecReRHyJORH99g+K8L1WSEORfzZ
i+2zYvCSMdtnxdT3jWyfFeH6rDCPYv4cyvZZMecd9yXuU9yf+WyfFRNPX7bP
iqDPqhvbZ0XQZ+XG9lkR9FndYfsGXfHf+5DtG2T2z1C2b5Dhm45s3yDh+gYR
1+N+82L7BgnXN4i4CfGSNts3yNR3Y7ZvkKlHDmzfIJM/t7N9g8x5V2P7Bpn9
GcD2DTLxzGb7BgnXN4j+TfRtcnywrv/wwZJ/+GDJP3yw5B8+WPIPHyz5hw+W
/MMHS/7hgyX/8MGSf/hgyT98sOQfPljyDx8s+YcPlqAPdghb9yaoe6uydW9X
1J8T2Lo3Qd3blK1LE9Sl29m6PfP7l9i6vSvq5NFs3Z6gbs/Ryck/dHKC35PR
OaWf8w8dnvmenL4A83wHtk5OUCePZevkBHXyUrYOT/6hwxPU4f3ZOrbrP3Rs
gjo2RycnqJNfZevk5B86PMHns/9+PfOcSrbO7/oPnZ/7fBmuDoz3w/+3biyz
5H/rzDIe/1uXlnH93zq2zJn/rXvLfPjfOrnMkP+tq8sU/W8dXmbP/9btZRT/
t84vs+B/9wVkKv93H0EG+7ScvgPzz6aXj8D2/f+nA6NP/TFHN0ZfNfqs0X+I
PuAXHF0afavdODo2+iwHcHRv9AVu5ujk6GOz5ejq6Lvy5ujw6BOaxNHt0dey
jqPzow9jAacvgL6BaZw+Ava5z3L6DtiXxThj/NlzFP/pwOj7H87RjdGnjr51
9KujrxrxLPp10Qf8nqNjo2/1E0f3Rp9lE0cnR1/gBY6ujj42wtHh0Xcl5uj2
6BOaxdH50dfSxukLoA/jCqePgL6B2Zy+A/a5Mc64/9lzKf/pwDhHsY+jG6Pv
P5KjM6NPvYWjS6OvGvkZ+qvRB4y+YPTDoG91NkcnR59lI0dXR19gHkeHRx+b
KUe3R9/VBY7Ojz4h9A2hXwh9Lb84fQT0Ybzh9B3QN4BxxvzDnvP5TwfGuRTE
m+h3xTkKxPvop0Lf/xSOLo0+dROOjo2+avRZo78LfcCVHJ0cfau3Obo6+iwH
cHR49AWWc3R79LH95Oj86LtazekLoE8IOH0E9LVYcPoO6MPAOGP+Z89N/acD
45zPVY5ujHMpfTg6M85RIH9FPyH6/lE/wHyCPvVsju6NvupAjk6OPuBIjq6O
vtUjHB0efZZJHN0efYHzOTo/+tjsOH0B9F1t4vQR0CfUwek7oK8F44zviz2H
9p8OjHNTyIcwn+CczxKOzoxzKe4cXRrnKK5wdGz0/d/g6N7oU9/B0cnRVz2c
o6ujD7iZo8OjbzWVo9ujzzKFo/OjL/A5py+APrYaTh8BfVd2nL4D+oQwzoh/
2HN9/+nAOIc2haMb49yUKUdnxjmfLRxdGudS5Fj6FqN7M3MVmM/R9084Ojn6
1G04ujr6qudydHj0AR/h6PboWzXn6Pzos2zm9AXQF3iI00dAH9seTt8BfVcY
Z8Sf7DnJ/3RgnOtbwtGNcQ4N9RKcR8C5KZyjwnyOcz7DOTo2zqXgnAriQ5yj
iOLo5Oj778nR1dGnfpOjw6Ov2oij26MPmMfR+dG3+oLTF0Cf5QdOHwF9gbac
vgP62DDOiP/Zc6f/6cA4J+nE0Y1xrm84R2fGOTQPji6Nc1PA0bFxzgd1b5z3
wbmUuRydHOco5nF0dfT9a3B0ePSpE45uj75q9Flj/UUfMPqCEf+gb/U0p4+A
PstYTt8BfYEYZ+Rf7Dne/3RgnDu15ujGOCeJc5M4L4lzfTs4ujTOoaGOjfNo
ODc1k6N745zPMo5OjnMpOKeC+QfnKFQ4Ojz6/odwdHv0qffl6Pzoq57M6Qug
D3gDp4+AvtVITt8BfZYYZzxf7Lno/3RgnOMdwNGNce60gKMz45ykCkeXxrm+
zRwdG+fQDnB0b5ybwjkqzD8452PF0dVxLkWZo8PjHAXOVSD+Qd//Wo7Ojz71
rZy+APqq33P6COgDzuL0HdC3inFG/YE9Z/6fDoxz0UEc3RjneG9ydGacO33H
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
snsx0b28/+21SRC9F0N2LzZT9F4M2b3YWdF7MWT3Ytmi92LI7sV6i96LIbsX
qxW9F0N2L3ZH9F4M2b3YKtF7XmT3vKqi97zI7nmn0Xtee/6eF9k972HRe15k
97yNove8yO55JUTveZHd8/YSvedFds/7WvT5B9nzz1LRe15k97xHRO95kd3z
jqf3vMjf8yK75w0VvedFds87V/SeF9k97yDRe15k97ymove8yO5574i6BWRu
wVzULSBzC+tE3QIyt2Aq6haQuYXfom4BmVvoS93CSd4t4CjqFvRE3QIyt2Ao
+jyP7Hl+qahbQOYWZETdAjK3cJm6hXe8W0DmFhRE3QIytzBA1C0gcwuyom4B
mVvYIuoWkLmFl6IOB5nDuSrqcJA5nN2iDgeZw7kn6nCQOZyTog4HmcMJFnU4
yBzOV1GHg8zhHBb9fIrs8+kVUYeDzOE8EHU4yBzOb1GHg8zh2Ig6HGQOR0PU
4SBzODtFHQ4yh7Nc1OEgczhNoq4MmSuzFHVlyFxZiagrQ+bKflBXtoF3ZShH
XVmJqCtD5socRF0ZMlemJOrKkLky0b3m/+21+btOxJUhc2WXRV0ZMlcmI+rK
kLkyVerK+vCuDOdTV7ZE1JUhc2W3RV0ZMlemKOrKkLmyL6JOEpmT9BN1ksic
5AhRJ4nMSUqLOklkTnKuqJNE5iS3izpJZE5yq6iTROYkRffE/7fXRlrUSSJz
kmKiThKZk0TqJNV5J4mx1El6ijpJZE7ST9RJInOSbqJOEpmTPCTqJJE5yVxR
94vM/aaIul9k7vchdb9fefeLitT9dqDudwfvfpG539ei7heZ++1F3e823v1i
CXW/+qLuF5n7jRY9D0d2Hq4g6n6Rud8S6n6P8+4XW6n7VRZ1v8jc729R94vM
/S6i7vcp735xDnW/ztT9WvLuF5n7taDutwfvftGBut+eoo4dmWPfTx37Z96x
4wbq2JE69iW8Y0fm2NOpY1fiHTsyx64r6thR6Ng3844dj1DHPoo69iO8Y0c9
6tjZ/c5cutdmLr3fERd17Mgc+35Rx47MsU+hjv0j79iROXZXUceOzLE/E3Xs
yBz7N+rYXXjHjieoY+9IHfsn3rGjE3XskqJzGcjmMvqKzmUgm8swpXMZsvxc
BrK5jCTRuQxkcxmX6FzGb34uA9lcRhydy7jPz2Ugm8v4LDqXgSvoXAbbF99d
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
sh1nl7/MWEvRRdoDGnD46Y6sr9GQm5F5neVTdrnat2S5lMF37z6c9SAaznXe
s4rlKcInbD3GchQaCYuuttGLgZ4bbANZbuhS7Lp2LC/Ur+xag6trYsB1Q0ma
NLfVXLk5y2ttGLnDSs0rBro3e9KD5eY8XkXdZHm5top3TFSyY6BHnUUhyy0O
vNbDjuUV13Q7yXKicPbpiacsN3qz7fkvLC+acMlUukO/YLGQ2xW+VkDJd3/p
bvpr+0ATlpvevVpdmpfWfrtIuvNtsWyeHcutl7lcXcby6kctm0l3qW8fShzD
egPz7S9FsL7AwZvzpTvKG+7OntGs+iGcrVP9yfoamoYrpbu/A01ep7HejMmY
jA2sL9P02yXpTu3ktHvHdFt6Q/KNgtWsD/UzJkO6q7o0bXNAWNUyY/clQn/q
+PvlX5VyIkD/1dMClvdxem0wj+V8HPM6L/2xMQLcdauesHzK331KlSyXsnjJ
d8ftwx9B+ut9QSyH9dki5gLLX001T9a55PQIxnqv0mQ5OL1Y/I/l30YUl3yL
KH0Ekzu+2s5yiJ0G96lk+cOBz0JuDDaOBPO+lV4sB+p+LWkvy38qDIlnuVuI
LHVSZTnctfk2a1n+1i68F8s9w95mQg5a+FoB1ao20t368vQbziyHrr+giTR/
PnqugnRnfdPqcx1ZD+C/VuttWP5/5MwnW9gueMcXRarNfkTDAlujNax/cc1s
ynS2Y511d0ca68E4r9Y4yfovXo+HFbPd5Y3N53Zc2fkuRK9vo836Rx8WjpHu
BJ+svmah1ssbAtqeUGE9Mrsv36S7tq5T5xveaLLcuGCs0Dub4xDqOycwGHrd
e+3N8lPRk0+NZ7mpORurFn5aFgzPBljGsLzPgl6H37Gcz+WhKuXNaoLh259r
df26FUBwre5mlmfLeNfJ6gA8gGfPn9mxXOE0FYPBLE9oU/Zo5bCDD2DsLj1v
lussfV6ozfKc10yuH/CIegBZE18bsFztjodRx1me9rzaRZZjhpTZFv1ZrnmY
z1hVlmde0tWe5cihadsO0ly58FXy+evtYDl+UPDY7cZy/avvXGnA8vyhIb+f
sp37S3VFv1ivQjvy6A3Wp/hze4sP22WPTrhZwnotU4u+vWZ9ltuwQ4/tiLtO
VItmvSKl8hZ7WJ/I6WzkS7Z73eugwZhfmnehWac9w1if69SihK1sp3lysIta
gp43GJiuu876dxHvekp3he/v/ZDdteNy43d65tK+3vri/0a0rvGDnYm2R1ke
zWrylTiWQ5voWmV67Jwf3BkWe4Xlp4oHKiDLTX2/YHbzywJ/6HFnkCfLCX4r
chzI8oGW/v6NDjz0h3cHI4NYTnPyH38Tls/0z/+1SaldAER86HKb5WQ/nBuo
yPKxNu7NjqRZB0DvR0pjWU4528G4hOWTE0YaLe8dHABLZ+z2YTlxBy/lT1LP
cIURy+WDeZ6Q0xe+Sr6fPTqzXgSMz218k/UkFp14MIL1IybkFX9z3xgIJeOy
AllPJWm1jTHrpzQefP3yuSuBcG7gEA3WEzKNOxTP+kF9z7ofZzvu1UGPbrCe
VtCQE0NYP2vixEBbtjs+b/mRPawnF2KtrMb6cYnHVnZnO9lThsWFtDHyhgdG
V4az3qKvmY1013lg6l63X72WG6sPMpD2HDuYdziCCm7wNLo8luX7lP20dViu
r3nVq3NNVrmBxygbTZavHBtzYCrLVSrfy7/247EbtFcc95HlW6tvZ89ludaQ
r3G5z7u5w5B7MZ4sX+xxq/silitequvk3HCLO1RM1Kxg+e67ak7hLNfdK3ps
zu0odxi3drg2y9cvaZgrzdX7zv5wUkvJAyIjhH6D8LUS0tT6e6aN94BzPY84
sX7J0fKJ0l5J8ufVubX/eUDw/F43WL/nw6Uhg1mv5+ejY5X9bnrA0Q2OY1m/
ynrJpL+sVzXm6bB9aU88wMGjy0HWb0u00bzLem3HOp5cx/bd9+++seL3qLuw
1uTPINYrdLg86ivbI/dde/Tkm7He4OfcWNrrnO/z1o7tZy/N87vZu/9y46xL
Qg/0acGtT7qTCiDkXccXxQvyQGfznTctfpUA3no20u5UAWwO8vdbeDMHLrTy
7ir1CXz6MCAvrQBcm2msHPIpC4ZMPGfAfA7xscpfo1oUgmXKAxMF9UzwmWr8
wDW5FDS+FwWljCuEbho23+zHpAM2zlJnPpPahm3XH99RCNo30tqYhD8Bq+Mr
BjKfzIrZjcYcSykEdUNU9i19DONyDwT+vlAGg2yMPqmWFcL8MaNSsSQKNken
HbgeWgbGCae2lCoWgem8J977w8KguPd+JeYNa/FobtOsTkVg1Xy3/+TjkvNt
QNuTph/LQP/vbnV3jSK469Vwt5mzN4zPy7n/RHLf1HXyle7izDDV8Zu9cxP0
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
EfuG/i4JhN/DTT8azHMFswl99jJPcvuWhv36RrvDgKWOi+fNyoDQ0N+jmZdP
a1y/3Wz3rtc81cJ2QwphUbi+PuPCHzqkSHfpRkQ5VR9aaGX8aYCXlBuOcii4
z/bzrKbs9mUccpDe3G+MP/7Z+/oR26vrqXN0E+OQsdbp5xh/dEwLlu4Xzo1b
9rNuwmb49FdJykPJz0K+FvK0EKckbkm8krwY5MkgPwZ5HMjrQD4H8g6Qh4D8
A9STp9489eWp1009b+p3Uw+ZesnUR6beLPVoqT9LPU/qfVLfk3qJ1FOkfiL1
6KhXR3066n1RD4z6X9RTot4S9ZXIX0M+G/LYkO+GPk/y3pBPhPwi5BUh/wX5
MMiDQb4G8jeQt4H8AuQbIM8A9eGpH0+9eOpvU5+betzUN6b+MfWOqR9LfVnq
yVKfk/qd1Ouk/iH1EamHSH056s9Rb476XdT3op6XfB+pEKiXRN4f8gCR/4c8
QfQ5ky+I/CzkayFPC/lEyC9CXhHyX5APgzwY5GsgfwN5G8gvQL4B8gxQH576
8dSLp/429bmpx019Y+ofU++Y+rHUl6WeLPU5qd9JvU7qH1IfkXqI1Jej/hz1
5uT7XYVAPS/yJZE/ibxJ5Feiz5k8S+S7If8NeW/Iz0K+FvK0kE+E/CLkFSH/
BfkwyINBvgbyN5C3gfwC5BsgzwD14akfT7146m9Tn5t63NQ3pv4x9Y6pH0t9
WerJUp+T+p3U66T+IfURqYco35crAurNkWeKvFPkmyIvFX3O5KcifxD5hMgj
RL4b8t+Q94b8LORrIU8L+UTIL0JeEfJfkA+DPBjkayB/A3kbyC9AvgHyDFAf
nvrx1Iun/jb1uanHTX1j6h9T75j6sdSXpZ4s9Tmp30m9Tvn+YRFQD5H8XOTr
Ik8X+bzocyavF/mYyM9EXibyB5FPiDxC5Lsh/w15b8jPQr4W8rSQT4T8IuQV
If8F+TDIg0G+BvI3kLeB/ALkGyDPAPXhqR9PvXjqb1Ofm3rc1Dem/jH1jqkf
S31Z6snK9zmLgHqd5DUjzxn5zciDRp8z+dDIb0W+K/JckY+J/EzkZSJ/EPmE
yCNEvhvy35D3hvws5GshTwv5RMgvQl4R8l+QD4M8GORrIH8DeRvIL0C+AfIM
UB+e+vHUi6f+NvW5qcdNfWPqH1PvWL4fWwTUk6X+obwHrlLmfSMPHPnfqP8p
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
OY9VuByMKybsyvf6rwh0lYaov4grk/lZyNdC3/9jFO0v39YLRONPKX+6ZKeA
kX1d320DKmH0pTenFUr9MGfMQ8N35kXgN9ellPVSk7hnhP8ckv38oZ4q9Vap
ryrv4ysE8vJRf5X6rNRjXcD7q8D7rNRjpZ0DOofT3gH1VKm3Sn3VUfLnTNkO
Qn2/Tqz8+RPrOX/K/PN/5c+fWM/5E+s5f8p81D7y508Unz8zhfOnzPe7Sf78
KfPTNpQ/f2I950+s5/wp81Xqyp8/ZX5F0f0I/3/3I2N+PzIW7kdYz/1I5vmv
4fejDcL9SOalD+X3o7PC/QjruR9hPfcjmac6Xv5+hPXcj7Ce+xHWcz/Ceu5H
sn352fL3I5nHUn63+5nMuzhC/r4v2wdJree+78Hv+9XCfR/rue9jPfd92S7A
Vvn7PtZz35d51/3l7/syT/g3+fs+iu/7KNz3sZ77PtZz38d67vtYz31fti8v
v4P+TOa9rOPvVy2F9yvZbsv3et6vjPj71Vjh/Uq2G9KSv1/tEt6vUPx+ZS28
X6H4/eqq8H4l2xFYJP9+hfW8X8k87Tby71dYz/sVit+v4oX3K5m3WVX+/Urm
Ga6Tf7+SeXH95N+vsJ73K5l3VPQeK9vTGfS/e4/Fet5jZTsj9+XfY2W7GDvl
32OxnvdYrOc9VubJ/yL/Hov1vMdiPe+xMm/2Yv4eGym8x8o8z8Pl32Mx5N/v
sTKP7nL591iZ91XEF/D/yBewHr4g23k5KM8XsB6+gPXwBayHL2A9fAHr4Qsy
D/w1eb6AYr5wQeALqM75gpk8X8AHIr5Q91nKF2Qe46PyfEHm3RXxR9n+VH2c
UcQTZbtFIm4o29npJs8HZbswIg4o2zEJled9st2NN/JcT7YTQfxulsDv8DTn
dyJOJ/Pwd5fncTJv/Fx57ibznJ/ifG2mwNdkXm5PeY4m80iLeJnMeyzivzIP
LfH3wQJ/l+2L1cd/6+PyIv4r86+KeL1sl0rEf2UeURHHl+0rifjv//gw5fm+
bCdIxH9lXkcR95ft3Yj4r8xPKMoDyHZbRPxX5tkT5QRk+yMi/ivzxanL5wdk
Oxoi/ivbdxDlCmR7EMry/Fe2c3pRPm8g2zUQ8V+Zh0qUQ5D5+aPk+a/Mp0T5
hAdCPkHmmRfxX5kXyEE+tyDzpYv4r8zjLcozyLzfonwO/h/zOSjO5/wn5HOw
nnzO/+xS8XzODSGfg/Xkc9D83/kc2U7NPp7PKRDyOSjO52wW8jkozufcF/I5
WE8+B+vJ52A9+RyZxx7k8zky73o0z5t9F/JmWF/eLFyUNzsq5M1ku3UonzdD
cd7stJA3w3ryZniB583M5fNmst0lI543QyFvhtE8b1YpnzdDcd5ssJA3w+c8
b1YgnzeT7YYslc+byXYu9sjnzVCcNzsl5M2wDc+bXZfPm8m899N4frKvkJ/E
/2N+EsX5yRAhP4mFPD+5lOcnNYX8JNaTn5TtiLXi+cmRQn5Stns1Rz4/+T87
TTw/+UfIT6I4P+km5CfxMM9PqsnnJ2W7Lbvk85Mozk+eFPKTOEuUn3ws5Cdl
Ow6z5fOTst0BygOPFfLAsj1Qo3rywOfl88BYTx4YxXng9kIeWLaLZyefB8Z6
8sBYTx5YtpN1l+eBJwh5YBTnge8JeWDZDtF9ngcuEPLAWE8eGMV5YBDywLhd
lAc+LeSBceq/88Cy3QcF+Xy7bKe14f8u3y7bzTwmn2+X7Twa83x7LyHfjvXk
21Gcb3cR8u1YT75dtlO2jefbXYV8O9aTb0dxvr1UyLfj3n/n27GefLtsF2aT
fL4d68m3y3Y3VvBeRoLQy5Dt54r6F7J9Ui/5/oVsT/OpfP8C6+lfyPYKQ3n/
QkvoX+Ad3r9oK9+/kO3BNZfvX2A9/QvZ3tZC3r9YKvQvcDbvX7jI9y9ke0b6
8v0LFPcvFgr9C9lejOgdVbZvUj6/Z1Htew883PPNuG8d72D8+Pm24Z8q4fXZ
4UoaqR7oquWUdb/nQ2x4s3OgYYbk/HzuRN+NNz1wa+XS54bHo7Hs0HSX1IBK
OFt7/uvYPR7YZcbppb2r4jB68O7maQcrwevGlGFdTT3Qc8XlI4XjkrHfpxKP
sGmVsHV/T+PjSh5oNibER+nAU3xSssmwS2vJr7N+3qG5W9xR58bbmVfNs/DC
lUHdo/ZWwDT7ecpdu7vjvm/BShamOag0d/Bmq/6S+7LvwSitODcMyn9gfV5y
L1528usnzYxy+DblQsthq9zQ7/MyrYiafOxkltT/6eZyWOncwfOLghv+KdW2
mhdWiAsv3bnAds+/yL8Dy3a6q+2i5zo9DsR1a8891/4RinNyB/p1TKqEcMtm
19d4B2Lfwfd6Nv8ahTd73/JR9q2EN4ct3RK2B6LGkp9PLiyMw32tMzra7K2E
H8ovM3uMCcQVs0o19txJwg6j5r+7NakSSs60TG+gGIjhnatfzq9MRTXjGQ5z
VSqhrmVWq4AeAXjrZBpYNMzBqy2GqetKztt3rDT2vk3wx6EOI34PLc5FNUcD
TwXJPcvv4MFFDdb4o3dEbEr30Hz0uHjguse2cnioHRQZoiT53/s59L2yqxCV
u3h+Y3vurYJUk0O8/HD1tvUuS7SeoXf2s4dsfzxt31SV6oAQ/OnUb5zagiis
nnFcZ51LJTR5pXFny+kQDJ+939Mg9jEe6fHIwJTttHbo/EfXJgSbzql5NaNT
EhYNylqlYFwJC7qMPXnCMAR37V5iVWGZiqZLnuj0bVQJbb5M3W9k/gAPWf9t
/uWK5OduyiR/w6flcMfHZm2KygOse9+w1tI2H3MdGp9Qkfx51Rj02/fncTD2
nhR4bsvQQnQd13DDJMnPvTID28dmu4JxTpu6B9+LizB87apbbCe95oR1caYF
YvtR3X3v9nmMGlqft9ktqITW23ranNRAnNJxsEfVrkS89njosNz+laBkf9Yj
oSQS31xfdWd+6BO0x7jlF6oqYMF2lY5LljxCLZVhJ9z65KNJtalal9WSzzNj
+LbQugiMXqdl0bGkAPOL1YvvNC2HJqo3HAdeicByT73UHm5F6H15hc6sa5Kf
t+rDs2tWxeKnwBfXGucl4N1mV8KnSL7P28xbFNd0XCy6Ol16OLrdE3w9dcb8
hPwK2DXfMTK/YwyqekW4JlwuwLD+6TfrasrAz2LJ78uF0XizUqdLlxlFmPus
NpftvNssTrkYoZSApsHrfk5Zm4IqS218SkMrYNwY2x29l8Zh8QrFjL5/CvG8
pWkx23/XOFllN29WBm5qXtKnb7Q7FuX3MWI9u+RtLd51npyBUzuMPTttWxAW
jG/eP+1RBWhbRPuGGGXgpWW3G5zUDsPWsc3nlt2ogCKPEWg+IAMLh22PDCuO
wrGzfrTqwbzcu1teVuuQgZcXtOz8oW8cJoU9+fTQqgKWqHcsm/A7HXsZdf5Z
bZGEPpsmnek5uAJ+6aRoqD1Lx72Nl05I35iKxpv3tEr+VQ4HJ+gveXs2HZ89
N1np55eJH0ccmV27vRysFawv7bNOx8jwLo+nn8tGqzYlene0y8F+6hlQHpWO
q9NPNatenYu7rKffKy+UnG8DPvV/ppyOuS9Mnat187GmX/yPWZLP4fS+IVYe
AWn4yHmCRt77AgyMfj/2P+0yeLfu4K2t49Ow/eY+utNPFmFvx3Y32N59dLdJ
bgMk99LqI32ybm0JwISiRHvWc3Rya7h7TYMsXNKq2NHpZQgG940M1reR3JtS
dhrrV2Ti2/vDbhuMiMKF3p59JkytAMWPJQ+jczOx76Hilg82P8ZTRzS2fuhb
Af2ywivHPczEa9rb9D+5JuKQbgueb68ph9/ubfY3ccnE42/jrqYEPcFDEaYR
Kx6XQ7hh98mlapk4saRl3teB2RiwVy3Y/k8ZmJ0+rVZZkoEef548yGyUix3/
+huPlPy9Ti26vL/8QQa+046ZY5uch8nbQj63WVAGjsaXd463y8BeyS7VMw8U
4MXV+3R8G0m+T4xfr5w9IQPtin9GdtEswgYGY+LP+kr+fTGobhwXko3eo4av
SD3/AO/OTZnDeqZtfq0Y0eBaNm6LrLvZtj3isk7R35o1q4Ap12dWnDshuTcZ
B+w5+CgWPSbMPnLnfTmYbtbxCtuQjQpFPee0epWAwTNvFG0IL4fnDxQhYko2
mmW9/K5TmYJf2oafuX2sHCw1ncd1dGK7OKlXf7vmoE+o4rC9e8rAo98wfZ8F
kv9+ou7aJ6s8NIs+pvRe8ufYpMXGIeN6ZuG1bx4BezoV4KWLI88r5ZXCwL8d
/ee/y8QfRuY/LLEQo5elNdqyrxSq1VMsms/MRcuGkWu2PXyEBiP1v7Oe71jb
Dxc/Ds/Fpw1KKz8rxWKDwS897CXn55md9F3/dMzFN7c/RrXvnYD5YUd9bHeW
Q4/8llFnqnMQM77ev6udgsvtj75tOLEcggszsmb9kJz7V33Ti/6ei3+8q/p0
/FkKGtpvH7mHZaPdsbBZkX756HDyWbPn10phKr4s77krG79m6G1ZaV6I4Q33
xxbPKIWvxa+t9CX3gVV3D3RyNYvBQ53H3GY963mha77vC85D0/6tu3mvjMe9
/d6rlQ4qh6n/Nf5RdVpyPjgV/rSJdTJ+f5ypnaZQDiMX93zycGIuGoXb4S7j
fJzUQX+G7aZS0OqZfby1Yi5Gz+xscfNdAfac9u7UwK6lMHPKvAFPTCWfY8GW
X9N941D5kfFi1nOf+2vtXuseBZj4yK9XhlcSRqwZvGdvWBlkKP2Y63FU8ucR
2n24/ZYCVLwxstM5yTl5zpRinXM+hag+YfUsv2eJ6DRjWQvmGSCe650k5bky
L/EV+dyazBtMvLjFECkvxvWcF3MvMe4QeYnr6+2Kcl8y7zTnmNBexDH57wfE
vx/OQ0FFxEPLBX4K0zg/zeX8lPNWUBXxVs5noaWIz3KeC9kinsv5L/wQ8V/O
i+GbiBdzvgxvRHyZ82hwE/Fozq9hsYhfc94NdiLezbk5jBNxc84xoerfHBPE
HHOkwENhKuehlBfl/BQWcX4azvkp563QRsRbqwQ+C1s4n13I+SznuZAs4rmc
/0KgiP9yXgy7RLyY82W4IOLLnEeDrYhHc34No0X8mvNuMBXxbs7NoWXJP7m5
zJ9N3JxzTJj+b44JYo45X+Ch0EfEQ50EfgpPRfy0p8BbYbWIt04R+CxEi/gs
z2/DLc5zNTnP5fwXeon4L+fF8FzEizlfBlsRX+Y8GhaLePRxgV+DkYhfc94N
v0W8m3NzaP1vbi77nImb3xI4Jrz5N8cEMcfkPBRui3hoW4GfwgzOT5U4P+W8
FRJFvFVX4LOQKuKzqwWeC41EPJfzX+gh4r+87wC/OS/uxXkx58twXcSXOY+G
tSIezfk1HBDxa867YYqId3NuDt3+zc1lnzNxc84x4d2/OSaIOSbvp8AVzkOp
n8L5KURxftqM89PfAm+F0yLeyvkszBLxWc5z4a2I53L+C+NF/JfzYngm4sWc
L8MREV/mPBrqRDya82v4IeLXnHeDlYh3c24OT/7NzWWfM3HzEwLHBIWqf3JM
EHPMWQIPBR3OQ6lvxfkp7RXL+OkNgbfCYxFvXSbwWdrrxt2cz64UeC604jw3
n/Nczn9hnYj/cl4Mp0W8mPNlCBPxZd53g3acR/fgPJrza2gr4tecd0OIiHdz
bg4N1/yTm8s+Z+LmXgLHhDWcY2bLc0wQc0zOQ2GWiIfOEfgpdOP81IXzU85b
oZrzVnPOW3cLfBZ6cT47g/NZznPBXMRzOf8FGxH/5bwY9ES8mPNlqBHxZc6j
YbaIR3N+De1E/JrzbtAS8W7OzcHz39xc9jkTN+c8Dgb/m8dB6L95HNTD40DM
4zgfpN0VGR9cL/BBEPNBzh9hyb/5I9SI+CPvt4IS55tN5fkmvBPxTb4XBU04
PzWR56egLOKnKgKfpV0vVJHns6Al4rO0s2bzb/4LPUT8l/NlGP1vvgwKC+T5
8i2BX8Mezq9vyfNrcBLxa87HYcm/+ThMF/Fxzt/h47/5O8wS8XfO92HXv/k+
pIj4Ps8bwPp/5w1gzL/zBrLvW1HeQPbfKW/QSuCYMPPfHBPEHJPzUFgo4qGc
n8JvET/lvBUeiHgr57OwVcRnuW8BijjPreU8l/NfaCziv5wXg5OIF3O+DLdE
fJnzaOgt4tGcX8MhEb/mvBtMRLybc3Oo/Dc3l33OxM05x4Qz/+aYIOaYuwUe
Crs5D33LeSjnp6Ap4qfc7wEPRbz1o8BnYYiIz3KeC4NFPJf7RiCG81/yjXBe
DBoiXsz5MtSJ+PIHgUeDk4hHc18KDBbxa8674bSId3NuDoov/snNZZ8zcXPO
MSH+3xwTxBzTSeChcFrEQzk/BUMRP+W8FVaKeCvns+Ah4rOc50KtiOdy/gvX
RfyX+3bAkPPisZwXc74Mh0V8mfNoGCDi0ZxfA4r4NefdkCfi3Zybw7B/c3PZ
50zcnHNMOPNvjglijjld4KFgwnko+ZE4PwUXET/lvBWKRLyV81kIFvFZ7neC
SZznAue5nP9CoIj/cl4M00S8mPNlMBHxZc6joZOIR3N+DZNE/Jrzbugg4t2c
m4PZv7m57HMmbs45JuRyjonyHBPEHJPzUFAV8VDOT0FbxE+HCbwVyjhvJZ+Y
q8Bn4Qfns/s4n+W+Mrgi4rmc/4K3iP9yXgxuIl7M+TJsEPFlzqNBQ8SjOb+G
7iJ+zXk3VIt4N+fm8O7f3Fz2ORM3byDwTSC+uYPzTc5DYZ2Ih3J+Cl1E/JT7
7sCK81by3VUKfBYacT57jvNZ7tODbZznOnCey/kv7BHxX86LYaCIF3O+DBYi
vsx5NKwW8WjOr2GDiF9z3g3eIt7N+TjkiPg4f/+UfZ70/rn37AXvyavcwG/i
VVfztZLv2yle/myPL2pB9CudODeoXLtu8lHJz4Evo/Uv1oaVg9sJq/xW3d1h
8JS73Q5E5aDJ/kEGbI/yUr63vv4Wd/hgtGJiQUYWDluuP4rtgU44PVytMsod
ft7eumxhVga+2qY4me2xYpMPDvZKHnBqW9WqLifS8JLjuGC2h6ulqJqxbI8H
vL3Q6KVO0wSs/mq5nu0+m98ct3nnTQ+IDOwyU+FaDN520jjMdrd7hN990CHV
AzwXfmkctzQCf2pNWch2z99dvTa+8L0HXAPdnzsX3sOmXQansN351EGPE3d9
94A+Nj6Hqxx8cOadYtfV3yrp/Rnuit6fpy99t+y4vz/YLG39oUVGHlpP23b3
iU85BOiZ2sQqBcCOOYPMXL/l4ORpRkvYXqeZZe9kq4UBsOfFKs2Dytn4e3yV
EttLtT/zn2/grQBof3hT76MqmTjky34Ltlc7ZmGQ0++qAGg6cWbANYV0nJV3
vAPbC+6k2zOr0jkQBrZ5fnSsdizGj+4QwfbH886+iAoPD4SZnQsd3LY+wpjH
nVPY/vubbpdU4nMD4du8TVEx2++j/rPZLV/nSO6tu6f3i/kQCMs3KXfbkHwd
e7kn7q36UAmtqtfq7a8IBPtTPTJjxh3EIqdhn8IqK8EsrN8I/40PoDq08tCj
Prl4tKjFDXhXDgcfFG8zD3wA30bePXzQMBsLLj2fwPZhB+w3P2xZ8gAaafj+
WTY2E6fmKr1SOV4Bbw/WGgzuEQLT/qimmhmko8GYj3ZsH3lOec2loQ9DoHrQ
pofT7CLR0tB8Adu19xjUwmppWghstxg9upVjMBY+T9OYjpVg3DxSrePzEMj5
u6CHjcJNfH1sw/CDmZXwzaC/o97bEFB7mGqsMt8R15vqrjpSVAm97h556KYW
Cf2M14XGrsnGn/69VrJ922v3zBfvXBMJa7t+mJW1NROXW4Vfm2ddAddzC5xd
7kZK7l9GhVUb09HDefhFtu/8Zf2ufUaeCDu3jm9ldekBPgporKnkVQnOhQOW
RN1DOLLg4Z7YQbfRaJ1/uynBldCyeH240SOEKrOg/VMGnsOzDyuO7pP8/h1n
fm/n8zFG8uupbHvslYlu2k9LtUHy92LHoLtP2sfCUxW1gNWe6dhYu80Qtjf9
YsbMz9PDYmHNp86a+839cfTolQWOZyrh1ItpP5pHxcL+Xn+fvHN3xkaHp+q9
dK2EO5X6xttj4yFsfr/f156ko2ra0Ndsz3pDq+rDS50SYOI1nz8pfy9j3Pnq
wOpNlbCW73q/Ee16t5+9b+H8v2lw7pOvWYpHPlrcyKku31AG79NXVqzWTYei
1f2V85/l4palQ3uxXddRb0drHViUDm0LdHO7NsvBrkbJJ9mu7sUb55q4HUuH
xn2tbwzXzMJloakd2K4xnO82a0VQOuwb5aJfop+B77M9Q9iutG3i8w4HstPh
mnKk6eKtT/FTu1M2bNdb+ZNJ6tbWGWAVve2IzfB4fHbHL5ft10963qd9eL8M
2GbjOuZETTTO7/+sXdn+Cni1vVX/kaMyYKC326gW88OxvWKDU5Wekp9LfYpe
RU7IgLgPL17F+N7FKMX7huphFTA1a8RUTcm97E5NrcK21KsY0D7jtV5iBYzT
dYhezPaMbvrtDwm3xXmOunlKkvtXfXvcVWMSXwzJyoDRaRHdzEpycU7f5xfZ
Hm58WaMB7oqZ8GTdlIqu6jmo0u9ZD7ZH/CZsmUX60ExoHqQ89qdxFobcuGnE
9qD3GR/wX7s0E6ZPexHqZZ6BPba07cT2uIcZ9y9UPZYJsxonLG394CnOamp4
nO2hbz9Q0en320w4b+Zw1WRBDK63mHJ57qgKKHFO2uD9MxN2j3o05OfgCOxu
pWDgukByPzqg1qxcKQs2XHE8+c75Hk5fn/ZswtYKCO2r1fRemyzoN8YvRTXr
GupovlDrJrkvlw1ufKRQcj9KCYp/Eep0ADWm/VD86iD5HPgu9li+iz2b7wIP
7nZE67R3Fth7BTxq2D8HV4XZGrBdZtvvvxI9c7JgjeHPswmzslArwWEY28X2
Oa4R+LJxNmwwmJTjvDYDQ4NvtbnQrRzm7l7++pNuNrxPcrh58/NTrI6HTLYL
P3Xd9mqF29mg57hKubrhIwzerLRtepsK+Krw10rykxGKze9kLdtzH8+/f3iZ
5YImX95XZZqQDYs2Df2mm3wDLYZklr6WfD6bx+9w3CG5jxTMDQ0eoOaAH21T
D/uMqah3b3rukZA4nzU5YP19sFXs0izsqeFzgO19t97o9trMLQdmT3h6/vX+
DCzb2dmX7a1D6FL/syk50G55F5spGmm4wuXKarZ3f2xrf9O7y3Ohp0nXrtlz
gzGmZSfr4qJy0L/q5We5SXJ+HjUu1O/OTbR45lrzU3Iv9lfq8eTq9lz442Y1
aZb+Wfw1OqzJ+W/l9LnTn4Ps87+08qeLrVYejFsS2MzjTAauuupbwHbhG1gZ
Tb8/Pw8M6v5E9JyVhrdmR51ohmWgrx3fZnOLfJhoPf2Kz7Hb2PpCaNzLy+Vg
4F64p0Q1H9rZV82+M+4C/ppdtN9Y8u9dfXvEK991yFJ9ng8aES13LNuRhmFW
Ted8lpzzF02f9ffL1QJYE2bSKWmICyrl9TjSfUa5bFd3rWhXl+evQJy/4n1b
2bmF+rac24KY277ie7Vv+F6tEt+r5fk3EOffOEcGNRFHXuqfpc08+GMy9M8x
L37J0tPKjFOHHdoxz+NoHoQ9bDfYfksBFC1s25lx7Qj1z7utexTA2TZV8zO8
kuBG98u7GAf/VvFCi3nhq2rmGM3wjYPO024sYNw8fEdbx9aKkvtgw8nTbr4r
gP61dmcYZ0/Kvse4PEyv3hC+yzgfknfUSbn8jNYdvlWdzgPTrzq1TayTocfw
w1qM46f3X/91X3Ae3PreYa73ynh4NMG+M+P+KvBtCfOwd1g+YRHzso9o9sOX
5QSqT8RU9tyVDVvznaxXmheCyYCqOJYreL1eHd3DskHH9v6kSL98yCxb24Ll
EN4lrsue9SMbPv/KGhj9PReGX/fvy3ILj/K8I89U54DeQ4uiu9op0O/L2dcs
57BNw9/lT0fJ9+3q5Kr2vRPA79hab5aL6N9/kfPH4bkwdOkh7SqlWPjo99iN
5Sj6nVkym3nVjb5d+ME866uvTPzKchenmx8MnP8uE553e/HBEguho3c/JZbT
UHw3ZNi4nlnwdeoV7z2dCmB3qaczy3Uo1iwa4bMgC5o4PHT7ZJUHRi1UlFkO
ZHKuzviOTlkwdOs5j9+uOfDoZag0N9L+YEvjiCnZMOB9lapuZQrouQ48zXIm
J71aeIZtyIbZf1z3tXqVAD3fPy9guRTfOd3Lz53IBu8f6QEHH8XC4n47DrMc
i0G2iX6Da9lQ2vt9b9X2CINjflWz3EuT8ppGcSHZsFk9v5p52QemTbBkOZl9
BXdWz5b8ezRa4/rtLppFMOTz6WSWq9H8k7J7vF0GLD5y6u3MAwXQ/kiqHsvh
jFYZdqD8QQY8jps2xTY5D4q021Wz3I6d8hX1ypIMCPE0CshslAvPz2masJzP
mU23WC4IXqsFZX4dmA1jdyZKc0GzTfbva+KSCe87mIekBD0Bkw6DwlmOaL//
14pxDzMh3X6KxSfXRGg48UsRyx3dOarEckrw2Wem/oPNj0F/19nNLKfUe66T
kX5FJkx6v/qFwYgoUClX0WC5Jr9Vw3ataZAFWxVVmp5/GQLrB829z3JQth3u
XmFeeLPmpd9vb5Gcl5s/OsxyUz8m3QzaOj4Nmt773HH6Scnfv6Mr/VjOKnPB
wFUeAWlg016ndd77AvjWaMBElsvaGJsz4JlyOsR/cbGv1s2HNirNfrMc183/
qk2UR6VDuaHH3y+rc0HR1iOY5b7G+01x2WedDpsN9cOnn8uG6k4LhrCcWKDp
LZYrg7an7i/188uE/cd1LFiuLDbEkeXQ4PX9kbPSN6ZCO/eSliyHduuga+mE
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
nb/PI3+fB3qfzxfe85G/5wO95/P3f+Tv/0Dv/2MFXoBKAi8A4gWHBL6AnC8A
8QXOI5DzCMn9ROARnF/gToFfAPELzjuQ8w4g3sG5CWoI3ASIm3DOgj4CZwHi
LJzLoLvAZYC4DOc4yDkOEMfpL3AffCNwHyDu80jgRHhS4ERAnGiiwJXwh8CV
gLiSs8Ch8L3AoYA4FOdWyLkVELeKFjgXVgicC4hz7Re4GHIuBsTFXvJ3pef8
nYnel5byPlUv3q+iXpXu2lK2d4ZJdhfPGsxzxXdr9KV7Z+7W15a3iQrA6qjd
112UM9CyIFSd7aNV8n5RDO8bUc9I7ULHzX6vQzDx0AmF1aHOOC3FUbq/5pqi
qW2RH4L6B9Vb3jnkj1HvFkr32tzzbns1aRaCAZG6DZQj03HZypMxbN+t9sO0
kaUZD/BeUsPtdjGZOK1VmR7bg/PkfR4F3u+hXs/55Zlsbw6nzYGvHzudw/Bf
y6R7c/PHb2D7dBh05paab9/bmPQwRbpPV/u657zXVxBLXJzv2jg9wPF7Pkj3
7MoHvU/peTMSVewTjlRtkPy79H5mHNu/+5LpVtliaSS2nXY+6O2WTDzRv6M/
28vr86VRm7GtI3HH4klpxTbZ2Pe06ka2rzeb93Oor0M9nZvmy9l+H1Zu6KV/
pdYBp7R8L93vU2/at3PCw1iMXj/MQDvRF40/PpXu/U3PjWD7gKi4c+LP0wuD
UVndWboP+MfzNtsTRMPY5QUVJpF41iBbuieYPrlXoJZ6LOYYq2R4q6fjn5Ty
4Wx/cG8P1/6+lTGodd8l7WrPTLy3t1kLtldYds5g2JfIGDR90Wh9g07Z6OdZ
0pntG55Qf/dxi0MMWmWv3x9fm4Pw1NOF7SFSn4f6PdTrmd3vL9tbxLmtLvf9
8PIA9ohoWMb2FrWPObB9Rmylf631yTXXUb/vLek+45hzTmzPEcMdpr1spHEf
o49pfmd7jmrnwtn+I37bfTnrR10E6mjGSPcf9yT+YnuR+Mtnzie3VzH4389t
0r3ImY22rFpSEI/N3RQ2vz6VhiqXt2uyfcmbvyxtdl6Lx48/9u9WfZKBTQ0d
PrM9ytlG8W0Pr4vHb6Pu+66MyUKlt67+bL+yvFXEhQa68XjmwkalTTdzULWv
3li2d5lgu9tzQ3kcDq48qLjBLg9bzZ72Z/bocnDh/SIz3jeinlGvuQPY/iau
36xTEx5ui7VLqqX7m02vrWB7ndj7VrcOiq19UC3OVbrXOXJSAdv3xPL2T//7
WnIXv5VZ57J9z3fhLYsydFPQ0XDgy00R4ag+vtNMtgc6TPUz2w/FyxNCzc1H
xeDliWM3sP1QjcKRUUdLk7Fg0u3sLkfiUeVLsXRvtO3sSWyfFK/qtPqz/OVT
/FuzWrpPatjnYrXOjGTst/OR9vEtGYimbs/Znmnbx8M6XW6ZjDXfJvxcY5OF
LRZl/mT7p+cmxtgGJSZhrUXuvQUzc/DPwPk6bC9VRSuoQ+N9SehxruPMs/3y
8O3ATvPZvuqoyk2R/XST8HN/71abK/Kx0zQvNbbHqsb7VNSvol4Vnd+En0+n
Iuo5v2E95zes5/yG9ZzfsJ7zG9ZzfsN6zm9Yz/lNlmtVlj+/yfw/ovObLGcb
JX9+k+V1M+XPb7KccK1at5bz22WhhYPpUvsul/FFzQAFtpM778t0tquLEz42
/Gg2Jx3rTgy3ZLu6CTmL2Q4vju+sd6Hu8AX85W29l+3w9t54Ky3wcTZqq29/
Yajghxbrouaz3V6/nqOXVvTOxrjZsbd1P6Zh47FXEtnO74UrYWwXGF062sdP
f5eBweNfSneBxy/IZDvCGGliN76s0VksNb0t3RGe8cJ3QRPJPTrixO4t2pdu
ot/UlDdsd3hXn8L4PpJ7dNio5g/HTw3Gd3mXLrOdYheVtIa/k3Lw/UP9TPXu
adh0+LgtbNf4v+Fv2A4yTtzRvbHt3gy0GWVexnaQG/pqHwlenYODSo6N2LM4
CzdavZfuJjeubTUqoF0+1to2vbVu6HH0s1ot3VmuzP3OdpnRdqlrr0H6N3DT
sn3N2C5z9J8Pxk0V8rFdktOQz03u43F3ozVsx/ntcfsjaRV5OOCq1tAgqwgs
2XWxHdt9PrMy3e7Dijx8fTQtt8/spxg4rKIZ24n+pXliaMjQPBxUMU9nUvsM
DLQuUWe70nWDf61OrMvFC38Sd/X4mokle+Y6sB3qxIY6jvbxuegMvQeGxWVj
oLuzBtut/mU/PPzatQJs9nJ4+/vn7bDA1Ee6c/1s0KN3szwLMNim8lnqK8nP
E0vLoWwXu79W9w+K5wvQ4YCSbmDcXfywckZftqPdNFU18qVdAaod6ZsU3jwc
51VYnmO72/EGx/arWBfg+qWpus3do9Hg/vwcttONLqNajfiaj49fnFv+62gq
em7PHcF2vXNLmra/FJOPKTFd2gzxS8f+J6JPsh1whZLKh5WO+Xi3etGAw5L7
dC/NmxVsN9xs2X7P2zPzsf3342dbzc3GgG/Bq9nOeKlJMdslx5pug5536ZKL
M4o+SHfJA/xvdnXXKEJDxXleM529cciMyGC2b76o+9BmWZ2K8G2IWu6k43cw
/rX/KbaH7pn8eVupYhE+bp74al9YGPa/46PM9tM736grVS0rxBGD/27Hkig0
q7E8zPbWI2Prxh1LKUSlgvGRN0of442dY+6yfXaFnu02Hd9RiBWW24sg/Ale
+qAxiO25L6zxv5cyrhB/3LN4aT8mHTudXtSd7b/7tDz5PapFIaYUndFVUM/E
mjubQ9levGZ0A928tAJUbxE6d8inLCyeam7E9uV1Nm02sDtVgKevRPsslPz7
Micppxvbo1+r5VuiO6kAX3utKyxekIdnbUa/Zfv1LrwXWsN7ojm8HxrLfejk
R6/kXmITwWMM3GOM5DHm3mPg3mMk7/Et7gfnnmQkT/JfwasM3KuM5FXmHmZY
LHiYkTzM3NsM3NuM5G12EDzPwD3PSJ5n7oWW+bnJC8090jBZ8EgjeaS5dxq4
dxrJO8091cA91Uieau61hk6C1xrJa8092MA92EgebO7Nhs2CNxvJm80928A9
20ie7feClxuaC15uJC/3VMHjDVaCxxvJ471K8H4D934jeb+HCZ5w4J5wDOSe
cO4VB+4VR/KKcw85VAseciQP+WLBWw7cW47kLeeec+CecyTPOfm4uRdd9vlz
jzpwjzqWcI96e8G7Dj8F7zqSd5172qGh4GlH8rSfFbzuwL3uSF537oEH7oFH
8sDfEbzx8F7wxiN547lnHrhnHskzz730wL30SF7644LHHrjHHsljz733wL33
SN577skHP+4Lp+//+4JXH34KXn0kr34TwcMP3MOP5OGfKXj7gXv7kbz9EYLn
H14Knn8kzz/fBQC+C4C0CzBQ2BEAviOAtCMwXdgdAL47gLQ7YCjsFADfKUDa
KeC7BsB3DZB2DeKFHQTgOwhIOwh8NwH4bgLSbgLfWQC+s4C0szCX+9DTuR+d
vOjUSxoh6iVRv0l8bqFekthbSP0j8bmIekZNRD0j6hOJz10tRb0hOqdR3wdE
5zp+DgTxOZD6O+JzI/V03onOmdTHEZ9LqXezUnSO5fsXwPcvkPYvFgt7GXBU
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
45rnR3LOv8YvP2lBzn2naeO69ibnvf/F98J5iNGREs7HtWjxXPid29mR6qFY
vQFLqD4Kjo87QXXRdrWKrs+JbvCG3Q+2Ex2RoTZ3ylGiH4YcDrgnR+7ddhtW
Oc0m93DFNzMDev/uqFV9bUTq7Etru1GlpO5KKfxIp/XWK345HiH1xT1mXX0c
qTfhd9dupXXGO+7UrX3kHPZQvFI3kZzLJ72GvKHnMdrklfD7XuV/W/i9T74s
/s6n+5oPofrs+sS8w4OJXuuu8dmM6rSSCo+KwUTHHKvK9vcguubnto5duxA9
s+2Nv+VkogOMCr4/GE90wagZeeeoHhiUqT6lDbl3bSateGpI7uFWdr9y6f0b
UNbW0ITcB7f3DBjYntwPrZqLY+m9EG3WmP6S1MGttpUTpUhd9FA3zaf1sFLn
nO06Ui/e130do0nqx/fy4taapG5o+Y90f0DO4fuAfsWzybk0U0ma9pycxw8V
0sJ5GHehmfqu8OObeC4cj3va/iD68uaBkFZUbxoW7EqmOrPbcskwD6LDju2Z
d2U60WX+R+/coXrs9pjXFtOJjmk1P31xGdE1KfZJWVTP2Fe7qGUR3bDJIkYl
jOiITirzYqh+SLh2tSe91xdVndvSi9zznbaZzKX3+811s6szyT13aeXHRXPJ
vTfF4OYAXXLfmRfsnbeN1PczfZNNXpJ67+Dbe0s3UucLQ023riJ18Gu3OaoL
SV2UOTVIQuvhesOJMWtIfXGQNfRbQepNxQvbQlpnHpQ7hK4g53zL1ambtMi5
hzsBh+l539JePJ87FqgJ53VBuXhOBzr/kk0jejqx/JcR1dcZ4ScjqK7WmZX1
ypDozuEPvE2oDg189GEJ1Z8/24W3Dia67Xv7n1ZNRMdVTDea+4Dotzk66ZNi
iE5q3eao622imzwnHBlcQfSS1D0nz/NEx7icVihOIbombOD541TPbBz0wr8L
0SWyu3/1ozplQcTNhVSfXDMxLFQn971J62SrOeT+19h4fHAtufc/r2semUnu
uWsn/8z/Se690B5xT+h956S8cfYXch9YHWvKn0Tuh9CSpBx6L3jcfl9/iNTZ
6pAR69VJ3Q2eGZ9P660k8M/9xaTeefa+tdST1L+CNVsqaN2D8s4+hqROpU2t
WBhI6pbfEWt/Wq/Mfv4Q6tFizU9CfTKYIdal45N2xuwi/YbGl81naf/hN3Xe
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
nfj7seWiQqx9fEtahdSBDSuMX9K/98q2IQl2CoV4TU3SYbhjCRQ5frtMv3el
6srN03UKcXHwxTPBx5PA4qP3LLxWDwWVgQFDlArRKkdX5YJyCsx7Orptx8P1
cLxjG81wqUK8r27bVjM8DTqN6NHp0sr/5fSx+/nvvTxJMnd9aqsC1FhdlxGb
UAhmZRXj6b//QviD+m4Z+YgdNxosm14MQ4pLvtPvK8ppg2xXz3wc2fP4x73R
JVD59dkA+rs929/k7mkH8v9bHjgRYpwCkoVXvo7KrIdmr0VjktbkY23PyCkG
DWmQGrZn67mA/+UzirmMT4HnM0p5KDVs1s7D8S8zTM3ti+HS2z1h9Pcz9kRH
91tt83DXnKz0Q5ISeKX7eh89R6uOv7+sHZ6LaoulJ+8cng631JamKVb9L2dT
fF/1KfB3VgeaRJnf2ZSDGlKX7L8ol8Iwi1O/6TlVZTmtTAf91T+cd3zC+EfO
HdUq1M0NMszAHXLHD0RphIPNmFV2daTff5el/MNzcgY+ryhNeb4iCozm7lu8
XUkCTx55lqqrZuDRhf1Mndsj3DBbtewZ0YHLh1oevN4pA9vMmtDLcUEyDOrt
OH/Cuwa41em98/h36ehS2+5MTtdUWGCVZvma3EeVcUV2vTEdTdt/tpPHdBh9
0fb1YI8GqGB5xOI71FnA36PuEtpOy354OtbtVrivFZ4LazMtcunnE7zS9uD8
hjQ8rWUYk7AmH0J2aB6h38vGjltfTgpPQ/vDvU6HShWC7A81Ffr7CXxk8eLS
pjTEgKNtXXyKoLuSWQ49Fxn9zOULVNOwcfwTSTvVEmi4MOU2rQ8De8r3nPHy
CZ7pZvGrMLAURnZUD6f1k/FveJNxoJx/G+WicthuXyp2+P7yWA/HKJDd81Qp
caEEDGLHtpZYp6JmyeILAzQRjKfe/xaoKYGmTzEzOhil4qdtC0zn2iaD26GN
WtvlJWCV2DOlZmwqukkFz64cmwpabc2rbBsaQP5L7x/5XVNxn87MsWPepYOZ
uvVG+nnyvGyml//q5O4bH5aGG6fg4TGDdU855EPHY5e6qZHf84YS6+W+yimo
7DV1q4pSIXTP/KFOz12P2Mj6Le8f44fuRoMPxRbBMumkQYWknoz7kRd0IPgx
vn3t7jJ/Tgk4W9WeoHV1ldVHD6VNj3Hc/LBxeqml8ErqdTS9dzif+pvxqpyT
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
37yOv/dl2uK9r/8Dlgup5A==
    "]],
  Axes->True,
  BoxRatios->{1, 1, 0.4},
  Method->{"RotationControl" -> "Globe"},
  PlotRange->{{-1, 1}, {-1, 1}, {0., 1.}},
  PlotRangePadding->{
    Scaled[0.02], 
    Scaled[0.02], 
    Scaled[0.02]}]], "Output", "PluginEmbeddedContent"]
}, Open  ]],

Cell["\<\
Ahora podemos realizar el mapeado del sistema (\[Xi],\[Eta]) al (x,y) \
utilizando la interpolaci\[OAcute]n de un elemeto de 4 nodos.\
\>", "Text", "PluginEmbeddedContent"],

Cell[CellGroupData[{

Cell[BoxData[{
 RowBox[{"x", "=", 
  RowBox[{"Simplify", "[", 
   RowBox[{
    RowBox[{"N1", " ", "x1"}], "+", 
    RowBox[{"N2", " ", "x2"}], "+", 
    RowBox[{"N3", " ", "x3"}], "+", 
    RowBox[{"N4", " ", "x4"}]}], "]"}]}], "\[IndentingNewLine]", 
 RowBox[{"y", "=", 
  RowBox[{"Simplify", "[", 
   RowBox[{
    RowBox[{"N1", " ", "y1"}], "+", 
    RowBox[{"N2", " ", "y2"}], "+", 
    RowBox[{"N3", " ", "y3"}], "+", 
    RowBox[{"N4", " ", "y4"}]}], "]"}]}]}], "Input", "PluginEmbeddedContent"],

Cell[BoxData[
 RowBox[{"2.5`", "\[VeryThinSpace]", "+", 
  RowBox[{"2.5`", " ", "\[Eta]"}]}]], "Output", "PluginEmbeddedContent",
 GeneratedCell->False,
 CellAutoOverwrite->False],

Cell[BoxData[
 RowBox[{"2.5`", "\[VeryThinSpace]", "+", 
  RowBox[{"2.5`", " ", "\[Xi]"}]}]], "Output", "PluginEmbeddedContent",
 GeneratedCell->False,
 CellAutoOverwrite->False],

Cell[BoxData[
 RowBox[{"2.5`", "\[VeryThinSpace]", "+", 
  RowBox[{"2.5`", " ", "\[Eta]"}]}]], "Output", "PluginEmbeddedContent",
 GeneratedCell->False,
 CellAutoOverwrite->False],

Cell[BoxData[
 RowBox[{"2.5`", "\[VeryThinSpace]", "+", 
  RowBox[{"2.5`", " ", "\[Xi]"}]}]], "Output", "PluginEmbeddedContent",
 GeneratedCell->False,
 CellAutoOverwrite->False],

Cell[BoxData[
 RowBox[{"2.5`", "\[VeryThinSpace]", "+", 
  RowBox[{"2.5`", " ", "\[Eta]"}]}]], "Output", "PluginEmbeddedContent",
 GeneratedCell->False,
 CellAutoOverwrite->False],

Cell[BoxData[
 RowBox[{"2.5`", "\[VeryThinSpace]", "+", 
  RowBox[{"2.5`", " ", "\[Xi]"}]}]], "Output", "PluginEmbeddedContent",
 GeneratedCell->False,
 CellAutoOverwrite->False],

Cell[BoxData[
 RowBox[{"2.5`", "\[VeryThinSpace]", "+", 
  RowBox[{"2.5`", " ", "\[Eta]"}]}]], "Output", "PluginEmbeddedContent",
 GeneratedCell->False,
 CellAutoOverwrite->False],

Cell[BoxData[
 RowBox[{"2.75`", "\[VeryThinSpace]", "+", 
  RowBox[{"\[Eta]", " ", 
   RowBox[{"(", 
    RowBox[{"0.25`", "\[VeryThinSpace]", "-", 
     RowBox[{"0.25`", " ", "\[Xi]"}]}], ")"}]}], "+", 
  RowBox[{"2.25`", " ", "\[Xi]"}]}]], "Output", "PluginEmbeddedContent",
 GeneratedCell->False,
 CellAutoOverwrite->False],

Cell[BoxData[
 RowBox[{"2.25`", "\[VeryThinSpace]", "+", 
  RowBox[{"\[Eta]", " ", 
   RowBox[{"(", 
    RowBox[{"0.25`", "\[VeryThinSpace]", "-", 
     RowBox[{"0.25`", " ", "\[Xi]"}]}], ")"}]}], "+", 
  RowBox[{"1.75`", " ", "\[Xi]"}]}]], "Output", "PluginEmbeddedContent",
 GeneratedCell->False,
 CellAutoOverwrite->False],

Cell[BoxData[
 RowBox[{"0.5`", "\[VeryThinSpace]", "+", 
  RowBox[{"0.5`", " ", "\[Eta]"}]}]], "Output", "PluginEmbeddedContent",
 GeneratedCell->False,
 CellAutoOverwrite->False]
}, Open  ]],

Cell["\<\
Ahora obtenemos el jacobiano de la transformaci\[OAcute]n:\
\>", "Text", "PluginEmbeddedContent"],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{
  RowBox[{"(", 
   RowBox[{"Jmat", "=", 
    RowBox[{"{", 
     RowBox[{
      RowBox[{"{", 
       RowBox[{
        RowBox[{"D", "[", 
         RowBox[{"x", ",", "\[Xi]"}], "]"}], ",", 
        RowBox[{"D", "[", 
         RowBox[{"x", ",", "\[Eta]"}], "]"}]}], "}"}], ",", 
      RowBox[{"{", 
       RowBox[{
        RowBox[{"D", "[", 
         RowBox[{"y", ",", "\[Xi]"}], "]"}], ",", 
        RowBox[{"D", "[", 
         RowBox[{"y", ",", "\[Eta]"}], "]"}]}], "}"}]}], "}"}]}], ")"}], "//",
   "MatrixForm"}]], "Input", "PluginEmbeddedContent"],

Cell[BoxData[
 TagBox[
  RowBox[{"(", "\[NoBreak]", GridBox[{
     {"2.5`", "0"},
     {"0", "2.5`"}
    },
    GridBoxAlignment->{
     "Columns" -> {{Center}}, "ColumnsIndexed" -> {}, "Rows" -> {{Baseline}}, 
      "RowsIndexed" -> {}},
    GridBoxSpacings->{"Columns" -> {
        Offset[0.27999999999999997`], {
         Offset[0.7]}, 
        Offset[0.27999999999999997`]}, "ColumnsIndexed" -> {}, "Rows" -> {
        Offset[0.2], {
         Offset[0.4]}, 
        Offset[0.2]}, "RowsIndexed" -> {}}], "\[NoBreak]", ")"}],
  Function[BoxForm`e$, 
   MatrixForm[BoxForm`e$]]]], "Output", "PluginEmbeddedContent",
 GeneratedCell->False,
 CellAutoOverwrite->False],

Cell[BoxData[
 TagBox[
  RowBox[{"(", "\[NoBreak]", GridBox[{
     {"2.5`", "0"},
     {"0", "2.5`"}
    },
    GridBoxAlignment->{
     "Columns" -> {{Center}}, "ColumnsIndexed" -> {}, "Rows" -> {{Baseline}}, 
      "RowsIndexed" -> {}},
    GridBoxSpacings->{"Columns" -> {
        Offset[0.27999999999999997`], {
         Offset[0.7]}, 
        Offset[0.27999999999999997`]}, "ColumnsIndexed" -> {}, "Rows" -> {
        Offset[0.2], {
         Offset[0.4]}, 
        Offset[0.2]}, "RowsIndexed" -> {}}], "\[NoBreak]", ")"}],
  Function[BoxForm`e$, 
   MatrixForm[BoxForm`e$]]]], "Output", "PluginEmbeddedContent",
 GeneratedCell->False,
 CellAutoOverwrite->False],

Cell[BoxData[
 TagBox[
  RowBox[{"(", "\[NoBreak]", GridBox[{
     {"2.5`", "0"},
     {"0", "2.5`"}
    },
    GridBoxAlignment->{
     "Columns" -> {{Center}}, "ColumnsIndexed" -> {}, "Rows" -> {{Baseline}}, 
      "RowsIndexed" -> {}},
    GridBoxSpacings->{"Columns" -> {
        Offset[0.27999999999999997`], {
         Offset[0.7]}, 
        Offset[0.27999999999999997`]}, "ColumnsIndexed" -> {}, "Rows" -> {
        Offset[0.2], {
         Offset[0.4]}, 
        Offset[0.2]}, "RowsIndexed" -> {}}], "\[NoBreak]", ")"}],
  Function[BoxForm`e$, 
   MatrixForm[BoxForm`e$]]]], "Output", "PluginEmbeddedContent",
 GeneratedCell->False,
 CellAutoOverwrite->False],

Cell[BoxData[
 TagBox[
  RowBox[{"(", "\[NoBreak]", GridBox[{
     {
      RowBox[{"2.25`", "\[VeryThinSpace]", "-", 
       RowBox[{"0.25`", " ", "\[Eta]"}]}], 
      RowBox[{"0.25`", "\[VeryThinSpace]", "-", 
       RowBox[{"0.25`", " ", "\[Xi]"}]}]},
     {"0", "2.5`"}
    },
    GridBoxAlignment->{
     "Columns" -> {{Center}}, "ColumnsIndexed" -> {}, "Rows" -> {{Baseline}}, 
      "RowsIndexed" -> {}},
    GridBoxSpacings->{"Columns" -> {
        Offset[0.27999999999999997`], {
         Offset[0.7]}, 
        Offset[0.27999999999999997`]}, "ColumnsIndexed" -> {}, "Rows" -> {
        Offset[0.2], {
         Offset[0.4]}, 
        Offset[0.2]}, "RowsIndexed" -> {}}], "\[NoBreak]", ")"}],
  Function[BoxForm`e$, 
   MatrixForm[BoxForm`e$]]]], "Output", "PluginEmbeddedContent",
 GeneratedCell->False,
 CellAutoOverwrite->False],

Cell[BoxData[
 TagBox[
  RowBox[{"(", "\[NoBreak]", GridBox[{
     {
      RowBox[{"1.75`", "\[VeryThinSpace]", "-", 
       RowBox[{"0.25`", " ", "\[Eta]"}]}], 
      RowBox[{"0.25`", "\[VeryThinSpace]", "-", 
       RowBox[{"0.25`", " ", "\[Xi]"}]}]},
     {"0", "0.5`"}
    },
    GridBoxAlignment->{
     "Columns" -> {{Center}}, "ColumnsIndexed" -> {}, "Rows" -> {{Baseline}}, 
      "RowsIndexed" -> {}},
    GridBoxSpacings->{"Columns" -> {
        Offset[0.27999999999999997`], {
         Offset[0.7]}, 
        Offset[0.27999999999999997`]}, "ColumnsIndexed" -> {}, "Rows" -> {
        Offset[0.2], {
         Offset[0.4]}, 
        Offset[0.2]}, "RowsIndexed" -> {}}], "\[NoBreak]", ")"}],
  Function[BoxForm`e$, 
   MatrixForm[BoxForm`e$]]]], "Output", "PluginEmbeddedContent",
 GeneratedCell->False,
 CellAutoOverwrite->False]
}, Open  ]],

Cell["Ahora calculamos el determinante del Jacobiano:", "Text", "PluginEmbeddedContent"],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{"jac", "=", 
  RowBox[{"Det", "[", "Jmat", "]"}]}]], "Input", "PluginEmbeddedContent"],

Cell[BoxData["6.25`"], "Output", "PluginEmbeddedContent",
 GeneratedCell->False,
 CellAutoOverwrite->False],

Cell[BoxData["6.25`"], "Output", "PluginEmbeddedContent",
 GeneratedCell->False,
 CellAutoOverwrite->False],

Cell[BoxData["6.25`"], "Output", "PluginEmbeddedContent",
 GeneratedCell->False,
 CellAutoOverwrite->False],

Cell[BoxData[
 RowBox[{"5.625`", "\[VeryThinSpace]", "-", 
  RowBox[{"0.625`", " ", "\[Eta]"}]}]], "Output", "PluginEmbeddedContent",
 GeneratedCell->False,
 CellAutoOverwrite->False],

Cell[BoxData[
 RowBox[{"0.875`", "\[VeryThinSpace]", "-", 
  RowBox[{"0.125`", " ", "\[Eta]"}]}]], "Output", "PluginEmbeddedContent",
 GeneratedCell->False,
 CellAutoOverwrite->False]
}, Open  ]],

Cell["\<\
Luego calculamos el gradiente de la funci\[OAcute]n\
\>", "Text", "PluginEmbeddedContent"],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{
  RowBox[{"(", 
   RowBox[{"Jinv", "=", 
    RowBox[{"Inverse", "[", "Jmat", "]"}]}], ")"}], "//", 
  "MatrixForm"}]], "Input", "PluginEmbeddedContent"],

Cell[BoxData[
 TagBox[
  RowBox[{"(", "\[NoBreak]", GridBox[{
     {"0.4`", "0.`"},
     {"0.`", "0.4`"}
    },
    GridBoxAlignment->{
     "Columns" -> {{Center}}, "ColumnsIndexed" -> {}, "Rows" -> {{Baseline}}, 
      "RowsIndexed" -> {}},
    GridBoxSpacings->{"Columns" -> {
        Offset[0.27999999999999997`], {
         Offset[0.7]}, 
        Offset[0.27999999999999997`]}, "ColumnsIndexed" -> {}, "Rows" -> {
        Offset[0.2], {
         Offset[0.4]}, 
        Offset[0.2]}, "RowsIndexed" -> {}}], "\[NoBreak]", ")"}],
  Function[BoxForm`e$, 
   MatrixForm[BoxForm`e$]]]], "Output", "PluginEmbeddedContent",
 GeneratedCell->False,
 CellAutoOverwrite->False],

Cell[BoxData[
 TagBox[
  RowBox[{"(", "\[NoBreak]", GridBox[{
     {"0.4`", "0.`"},
     {"0.`", "0.4`"}
    },
    GridBoxAlignment->{
     "Columns" -> {{Center}}, "ColumnsIndexed" -> {}, "Rows" -> {{Baseline}}, 
      "RowsIndexed" -> {}},
    GridBoxSpacings->{"Columns" -> {
        Offset[0.27999999999999997`], {
         Offset[0.7]}, 
        Offset[0.27999999999999997`]}, "ColumnsIndexed" -> {}, "Rows" -> {
        Offset[0.2], {
         Offset[0.4]}, 
        Offset[0.2]}, "RowsIndexed" -> {}}], "\[NoBreak]", ")"}],
  Function[BoxForm`e$, 
   MatrixForm[BoxForm`e$]]]], "Output", "PluginEmbeddedContent",
 GeneratedCell->False,
 CellAutoOverwrite->False],

Cell[BoxData[
 TagBox[
  RowBox[{"(", "\[NoBreak]", GridBox[{
     {"0.4`", "0.`"},
     {"0.`", "0.4`"}
    },
    GridBoxAlignment->{
     "Columns" -> {{Center}}, "ColumnsIndexed" -> {}, "Rows" -> {{Baseline}}, 
      "RowsIndexed" -> {}},
    GridBoxSpacings->{"Columns" -> {
        Offset[0.27999999999999997`], {
         Offset[0.7]}, 
        Offset[0.27999999999999997`]}, "ColumnsIndexed" -> {}, "Rows" -> {
        Offset[0.2], {
         Offset[0.4]}, 
        Offset[0.2]}, "RowsIndexed" -> {}}], "\[NoBreak]", ")"}],
  Function[BoxForm`e$, 
   MatrixForm[BoxForm`e$]]]], "Output", "PluginEmbeddedContent",
 GeneratedCell->False,
 CellAutoOverwrite->False],

Cell[BoxData[
 TagBox[
  RowBox[{"(", "\[NoBreak]", GridBox[{
     {
      FractionBox["2.5`", 
       RowBox[{"5.625`", "\[VeryThinSpace]", "-", 
        RowBox[{"0.625`", " ", "\[Eta]"}]}]], 
      FractionBox[
       RowBox[{
        RowBox[{"-", "0.25`"}], "+", 
        RowBox[{"0.25`", " ", "\[Xi]"}]}], 
       RowBox[{"5.625`", "\[VeryThinSpace]", "-", 
        RowBox[{"0.625`", " ", "\[Eta]"}]}]]},
     {"0", 
      FractionBox[
       RowBox[{"2.25`", "\[VeryThinSpace]", "-", 
        RowBox[{"0.25`", " ", "\[Eta]"}]}], 
       RowBox[{"5.625`", "\[VeryThinSpace]", "-", 
        RowBox[{"0.625`", " ", "\[Eta]"}]}]]}
    },
    GridBoxAlignment->{
     "Columns" -> {{Center}}, "ColumnsIndexed" -> {}, "Rows" -> {{Baseline}}, 
      "RowsIndexed" -> {}},
    GridBoxSpacings->{"Columns" -> {
        Offset[0.27999999999999997`], {
         Offset[0.7]}, 
        Offset[0.27999999999999997`]}, "ColumnsIndexed" -> {}, "Rows" -> {
        Offset[0.2], {
         Offset[0.4]}, 
        Offset[0.2]}, "RowsIndexed" -> {}}], "\[NoBreak]", ")"}],
  Function[BoxForm`e$, 
   MatrixForm[BoxForm`e$]]]], "Output", "PluginEmbeddedContent",
 GeneratedCell->False,
 CellAutoOverwrite->False],

Cell[BoxData[
 TagBox[
  RowBox[{"(", "\[NoBreak]", GridBox[{
     {
      FractionBox["0.5`", 
       RowBox[{"0.875`", "\[VeryThinSpace]", "-", 
        RowBox[{"0.125`", " ", "\[Eta]"}]}]], 
      FractionBox[
       RowBox[{
        RowBox[{"-", "0.25`"}], "+", 
        RowBox[{"0.25`", " ", "\[Xi]"}]}], 
       RowBox[{"0.875`", "\[VeryThinSpace]", "-", 
        RowBox[{"0.125`", " ", "\[Eta]"}]}]]},
     {"0", 
      FractionBox[
       RowBox[{"1.75`", "\[VeryThinSpace]", "-", 
        RowBox[{"0.25`", " ", "\[Eta]"}]}], 
       RowBox[{"0.875`", "\[VeryThinSpace]", "-", 
        RowBox[{"0.125`", " ", "\[Eta]"}]}]]}
    },
    GridBoxAlignment->{
     "Columns" -> {{Center}}, "ColumnsIndexed" -> {}, "Rows" -> {{Baseline}}, 
      "RowsIndexed" -> {}},
    GridBoxSpacings->{"Columns" -> {
        Offset[0.27999999999999997`], {
         Offset[0.7]}, 
        Offset[0.27999999999999997`]}, "ColumnsIndexed" -> {}, "Rows" -> {
        Offset[0.2], {
         Offset[0.4]}, 
        Offset[0.2]}, "RowsIndexed" -> {}}], "\[NoBreak]", ")"}],
  Function[BoxForm`e$, 
   MatrixForm[BoxForm`e$]]]], "Output", "PluginEmbeddedContent",
 GeneratedCell->False,
 CellAutoOverwrite->False]
}, Open  ]],

Cell["Por lo tanto ahora tenemos las derivadas inversas:", "Text", "PluginEmbeddedContent"],

Cell[CellGroupData[{

Cell[BoxData[{
 RowBox[{"d\[Xi]dx", "=", 
  RowBox[{"Jinv", "[", 
   RowBox[{"[", 
    RowBox[{"1", ",", "1"}], "]"}], "]"}]}], "\[IndentingNewLine]", 
 RowBox[{"d\[Xi]dy", "=", 
  RowBox[{"Jinv", "[", 
   RowBox[{"[", 
    RowBox[{"1", ",", "2"}], "]"}], "]"}]}], "\[IndentingNewLine]", 
 RowBox[{"d\[Eta]dx", "=", 
  RowBox[{"Jinv", "[", 
   RowBox[{"[", 
    RowBox[{"2", ",", "1"}], "]"}], "]"}]}], "\[IndentingNewLine]", 
 RowBox[{"d\[Eta]dy", "=", 
  RowBox[{"Jinv", "[", 
   RowBox[{"[", 
    RowBox[{"2", ",", "2"}], "]"}], "]"}]}]}], "Input", \
"PluginEmbeddedContent"],

Cell[BoxData["0.4`"], "Output", "PluginEmbeddedContent",
 GeneratedCell->False,
 CellAutoOverwrite->False],

Cell[BoxData["0.`"], "Output", "PluginEmbeddedContent",
 GeneratedCell->False,
 CellAutoOverwrite->False],

Cell[BoxData["0.`"], "Output", "PluginEmbeddedContent",
 GeneratedCell->False,
 CellAutoOverwrite->False],

Cell[BoxData["0.4`"], "Output", "PluginEmbeddedContent",
 GeneratedCell->False,
 CellAutoOverwrite->False],

Cell[BoxData["0.4`"], "Output", "PluginEmbeddedContent",
 GeneratedCell->False,
 CellAutoOverwrite->False],

Cell[BoxData["0.`"], "Output", "PluginEmbeddedContent",
 GeneratedCell->False,
 CellAutoOverwrite->False],

Cell[BoxData["0.`"], "Output", "PluginEmbeddedContent",
 GeneratedCell->False,
 CellAutoOverwrite->False],

Cell[BoxData["0.4`"], "Output", "PluginEmbeddedContent",
 GeneratedCell->False,
 CellAutoOverwrite->False],

Cell[BoxData["0.4`"], "Output", "PluginEmbeddedContent",
 GeneratedCell->False,
 CellAutoOverwrite->False],

Cell[BoxData["0.`"], "Output", "PluginEmbeddedContent",
 GeneratedCell->False,
 CellAutoOverwrite->False],

Cell[BoxData["0.`"], "Output", "PluginEmbeddedContent",
 GeneratedCell->False,
 CellAutoOverwrite->False],

Cell[BoxData["0.4`"], "Output", "PluginEmbeddedContent",
 GeneratedCell->False,
 CellAutoOverwrite->False],

Cell[BoxData[
 FractionBox[
  RowBox[{"2.25`", "\[VeryThinSpace]", "-", 
   RowBox[{"0.25`", " ", "\[Eta]"}]}], 
  RowBox[{"5.625`", "\[VeryThinSpace]", "-", 
   RowBox[{"0.625`", " ", "\[Eta]"}]}]]], "Output", "PluginEmbeddedContent",
 GeneratedCell->False,
 CellAutoOverwrite->False],

Cell[BoxData["0"], "Output", "PluginEmbeddedContent",
 GeneratedCell->False,
 CellAutoOverwrite->False],

Cell[BoxData[
 FractionBox[
  RowBox[{
   RowBox[{"-", "0.25`"}], "+", 
   RowBox[{"0.25`", " ", "\[Xi]"}]}], 
  RowBox[{"5.625`", "\[VeryThinSpace]", "-", 
   RowBox[{"0.625`", " ", "\[Eta]"}]}]]], "Output", "PluginEmbeddedContent",
 GeneratedCell->False,
 CellAutoOverwrite->False],

Cell[BoxData[
 FractionBox["2.5`", 
  RowBox[{"5.625`", "\[VeryThinSpace]", "-", 
   RowBox[{"0.625`", " ", "\[Eta]"}]}]]], "Output", "PluginEmbeddedContent",
 GeneratedCell->False,
 CellAutoOverwrite->False],

Cell[BoxData[
 FractionBox["0.5`", 
  RowBox[{"0.875`", "\[VeryThinSpace]", "-", 
   RowBox[{"0.125`", " ", "\[Eta]"}]}]]], "Output", "PluginEmbeddedContent",
 GeneratedCell->False,
 CellAutoOverwrite->False],

Cell[BoxData[
 FractionBox[
  RowBox[{
   RowBox[{"-", "0.25`"}], "+", 
   RowBox[{"0.25`", " ", "\[Xi]"}]}], 
  RowBox[{"0.875`", "\[VeryThinSpace]", "-", 
   RowBox[{"0.125`", " ", "\[Eta]"}]}]]], "Output", "PluginEmbeddedContent",
 GeneratedCell->False,
 CellAutoOverwrite->False],

Cell[BoxData["0"], "Output", "PluginEmbeddedContent",
 GeneratedCell->False,
 CellAutoOverwrite->False],

Cell[BoxData[
 FractionBox[
  RowBox[{"1.75`", "\[VeryThinSpace]", "-", 
   RowBox[{"0.25`", " ", "\[Eta]"}]}], 
  RowBox[{"0.875`", "\[VeryThinSpace]", "-", 
   RowBox[{"0.125`", " ", "\[Eta]"}]}]]], "Output", "PluginEmbeddedContent",
 GeneratedCell->False,
 CellAutoOverwrite->False]
}, Open  ]],

Cell["\<\
Ahora podemos utilizar la regla de la cadena para obtener el gradiente de las \
funciones de forma:\
\>", "Text", "PluginEmbeddedContent"],

Cell[CellGroupData[{

Cell[BoxData[{
 RowBox[{"dN1dx", "=", 
  RowBox[{
   RowBox[{
    RowBox[{"D", "[", 
     RowBox[{"N1", ",", "\[Xi]"}], "]"}], "d\[Xi]dx"}], "+", 
   RowBox[{
    RowBox[{"D", "[", 
     RowBox[{"N1", ",", "\[Eta]"}], "]"}], 
    "d\[Eta]dx"}]}]}], "\[IndentingNewLine]", 
 RowBox[{"dN1dy", "=", 
  RowBox[{
   RowBox[{
    RowBox[{"D", "[", 
     RowBox[{"N1", ",", "\[Xi]"}], "]"}], "d\[Xi]dy"}], "+", 
   RowBox[{
    RowBox[{"D", "[", 
     RowBox[{"N1", ",", "\[Eta]"}], "]"}], "d\[Eta]dy"}]}]}]}], "Input", \
"PluginEmbeddedContent"],

Cell[BoxData[
 RowBox[{"0.`", "\[VeryThinSpace]", "+", 
  RowBox[{"0.1`", " ", 
   RowBox[{"(", 
    RowBox[{
     RowBox[{"-", "1"}], "+", "\[Xi]"}], ")"}]}]}]], "Output", \
"PluginEmbeddedContent",
 GeneratedCell->False,
 CellAutoOverwrite->False],

Cell[BoxData[
 RowBox[{"0.`", "\[VeryThinSpace]", "+", 
  RowBox[{"0.1`", " ", 
   RowBox[{"(", 
    RowBox[{
     RowBox[{"-", "1"}], "+", "\[Eta]"}], ")"}]}]}]], "Output", \
"PluginEmbeddedContent",
 GeneratedCell->False,
 CellAutoOverwrite->False],

Cell[BoxData[
 RowBox[{"0.`", "\[VeryThinSpace]", "+", 
  RowBox[{"0.1`", " ", 
   RowBox[{"(", 
    RowBox[{
     RowBox[{"-", "1"}], "+", "\[Xi]"}], ")"}]}]}]], "Output", \
"PluginEmbeddedContent",
 GeneratedCell->False,
 CellAutoOverwrite->False],

Cell[BoxData[
 RowBox[{"0.`", "\[VeryThinSpace]", "+", 
  RowBox[{"0.1`", " ", 
   RowBox[{"(", 
    RowBox[{
     RowBox[{"-", "1"}], "+", "\[Eta]"}], ")"}]}]}]], "Output", \
"PluginEmbeddedContent",
 GeneratedCell->False,
 CellAutoOverwrite->False],

Cell[BoxData[
 RowBox[{"0.`", "\[VeryThinSpace]", "+", 
  RowBox[{"0.1`", " ", 
   RowBox[{"(", 
    RowBox[{
     RowBox[{"-", "1"}], "+", "\[Xi]"}], ")"}]}]}]], "Output", \
"PluginEmbeddedContent",
 GeneratedCell->False,
 CellAutoOverwrite->False],

Cell[BoxData[
 RowBox[{"0.`", "\[VeryThinSpace]", "+", 
  RowBox[{"0.1`", " ", 
   RowBox[{"(", 
    RowBox[{
     RowBox[{"-", "1"}], "+", "\[Eta]"}], ")"}]}]}]], "Output", \
"PluginEmbeddedContent",
 GeneratedCell->False,
 CellAutoOverwrite->False],

Cell[BoxData[
 RowBox[{
  FractionBox[
   RowBox[{
    RowBox[{"(", 
     RowBox[{
      RowBox[{"-", "1"}], "+", "\[Eta]"}], ")"}], " ", 
    RowBox[{"(", 
     RowBox[{
      RowBox[{"-", "0.25`"}], "+", 
      RowBox[{"0.25`", " ", "\[Xi]"}]}], ")"}]}], 
   RowBox[{"4", " ", 
    RowBox[{"(", 
     RowBox[{"5.625`", "\[VeryThinSpace]", "-", 
      RowBox[{"0.625`", " ", "\[Eta]"}]}], ")"}]}]], "+", 
  FractionBox[
   RowBox[{
    RowBox[{"(", 
     RowBox[{"2.25`", "\[VeryThinSpace]", "-", 
      RowBox[{"0.25`", " ", "\[Eta]"}]}], ")"}], " ", 
    RowBox[{"(", 
     RowBox[{
      RowBox[{"-", "1"}], "+", "\[Xi]"}], ")"}]}], 
   RowBox[{"4", " ", 
    RowBox[{"(", 
     RowBox[{"5.625`", "\[VeryThinSpace]", "-", 
      RowBox[{"0.625`", " ", "\[Eta]"}]}], ")"}]}]]}]], "Output", \
"PluginEmbeddedContent",
 GeneratedCell->False,
 CellAutoOverwrite->False],

Cell[BoxData[
 FractionBox[
  RowBox[{"0.625`", " ", 
   RowBox[{"(", 
    RowBox[{
     RowBox[{"-", "1"}], "+", "\[Eta]"}], ")"}]}], 
  RowBox[{"5.625`", "\[VeryThinSpace]", "-", 
   RowBox[{"0.625`", " ", "\[Eta]"}]}]]], "Output", "PluginEmbeddedContent",
 GeneratedCell->False,
 CellAutoOverwrite->False],

Cell[BoxData[
 FractionBox[
  RowBox[{"0.125`", " ", 
   RowBox[{"(", 
    RowBox[{
     RowBox[{"-", "1"}], "+", "\[Eta]"}], ")"}]}], 
  RowBox[{"0.875`", "\[VeryThinSpace]", "-", 
   RowBox[{"0.125`", " ", "\[Eta]"}]}]]], "Output", "PluginEmbeddedContent",
 GeneratedCell->False,
 CellAutoOverwrite->False],

Cell[BoxData[
 RowBox[{
  FractionBox[
   RowBox[{
    RowBox[{"(", 
     RowBox[{
      RowBox[{"-", "1"}], "+", "\[Eta]"}], ")"}], " ", 
    RowBox[{"(", 
     RowBox[{
      RowBox[{"-", "0.25`"}], "+", 
      RowBox[{"0.25`", " ", "\[Xi]"}]}], ")"}]}], 
   RowBox[{"4", " ", 
    RowBox[{"(", 
     RowBox[{"0.875`", "\[VeryThinSpace]", "-", 
      RowBox[{"0.125`", " ", "\[Eta]"}]}], ")"}]}]], "+", 
  FractionBox[
   RowBox[{
    RowBox[{"(", 
     RowBox[{"1.75`", "\[VeryThinSpace]", "-", 
      RowBox[{"0.25`", " ", "\[Eta]"}]}], ")"}], " ", 
    RowBox[{"(", 
     RowBox[{
      RowBox[{"-", "1"}], "+", "\[Xi]"}], ")"}]}], 
   RowBox[{"4", " ", 
    RowBox[{"(", 
     RowBox[{"0.875`", "\[VeryThinSpace]", "-", 
      RowBox[{"0.125`", " ", "\[Eta]"}]}], ")"}]}]]}]], "Output", \
"PluginEmbeddedContent",
 GeneratedCell->False,
 CellAutoOverwrite->False]
}, Open  ]],

Cell[CellGroupData[{

Cell[BoxData[{
 RowBox[{"dN2dx", "=", 
  RowBox[{
   RowBox[{
    RowBox[{"D", "[", 
     RowBox[{"N2", ",", "\[Xi]"}], "]"}], "d\[Xi]dx"}], "+", 
   RowBox[{
    RowBox[{"D", "[", 
     RowBox[{"N2", ",", "\[Eta]"}], "]"}], 
    "d\[Eta]dx"}]}]}], "\[IndentingNewLine]", 
 RowBox[{"dN2dy", "=", 
  RowBox[{
   RowBox[{
    RowBox[{"D", "[", 
     RowBox[{"N2", ",", "\[Xi]"}], "]"}], "d\[Xi]dy"}], "+", 
   RowBox[{
    RowBox[{"D", "[", 
     RowBox[{"N2", ",", "\[Eta]"}], "]"}], "d\[Eta]dy"}]}]}]}], "Input", \
"PluginEmbeddedContent"],

Cell[BoxData[
 RowBox[{"0.`", "\[VeryThinSpace]", "+", 
  RowBox[{"0.1`", " ", 
   RowBox[{"(", 
    RowBox[{
     RowBox[{"-", "1"}], "-", "\[Xi]"}], ")"}]}]}]], "Output", \
"PluginEmbeddedContent",
 GeneratedCell->False,
 CellAutoOverwrite->False],

Cell[BoxData[
 RowBox[{"0.`", "\[VeryThinSpace]", "+", 
  RowBox[{"0.1`", " ", 
   RowBox[{"(", 
    RowBox[{"1", "-", "\[Eta]"}], ")"}]}]}]], "Output", \
"PluginEmbeddedContent",
 GeneratedCell->False,
 CellAutoOverwrite->False],

Cell[BoxData[
 RowBox[{"0.`", "\[VeryThinSpace]", "+", 
  RowBox[{"0.1`", " ", 
   RowBox[{"(", 
    RowBox[{
     RowBox[{"-", "1"}], "-", "\[Xi]"}], ")"}]}]}]], "Output", \
"PluginEmbeddedContent",
 GeneratedCell->False,
 CellAutoOverwrite->False],

Cell[BoxData[
 RowBox[{"0.`", "\[VeryThinSpace]", "+", 
  RowBox[{"0.1`", " ", 
   RowBox[{"(", 
    RowBox[{"1", "-", "\[Eta]"}], ")"}]}]}]], "Output", \
"PluginEmbeddedContent",
 GeneratedCell->False,
 CellAutoOverwrite->False],

Cell[BoxData[
 RowBox[{"0.`", "\[VeryThinSpace]", "+", 
  RowBox[{"0.1`", " ", 
   RowBox[{"(", 
    RowBox[{
     RowBox[{"-", "1"}], "-", "\[Xi]"}], ")"}]}]}]], "Output", \
"PluginEmbeddedContent",
 GeneratedCell->False,
 CellAutoOverwrite->False],

Cell[BoxData[
 RowBox[{"0.`", "\[VeryThinSpace]", "+", 
  RowBox[{"0.1`", " ", 
   RowBox[{"(", 
    RowBox[{"1", "-", "\[Eta]"}], ")"}]}]}]], "Output", \
"PluginEmbeddedContent",
 GeneratedCell->False,
 CellAutoOverwrite->False],

Cell[BoxData[
 RowBox[{
  FractionBox[
   RowBox[{
    RowBox[{"(", 
     RowBox[{"2.25`", "\[VeryThinSpace]", "-", 
      RowBox[{"0.25`", " ", "\[Eta]"}]}], ")"}], " ", 
    RowBox[{"(", 
     RowBox[{
      RowBox[{"-", "1"}], "-", "\[Xi]"}], ")"}]}], 
   RowBox[{"4", " ", 
    RowBox[{"(", 
     RowBox[{"5.625`", "\[VeryThinSpace]", "-", 
      RowBox[{"0.625`", " ", "\[Eta]"}]}], ")"}]}]], "+", 
  FractionBox[
   RowBox[{
    RowBox[{"(", 
     RowBox[{"1", "-", "\[Eta]"}], ")"}], " ", 
    RowBox[{"(", 
     RowBox[{
      RowBox[{"-", "0.25`"}], "+", 
      RowBox[{"0.25`", " ", "\[Xi]"}]}], ")"}]}], 
   RowBox[{"4", " ", 
    RowBox[{"(", 
     RowBox[{"5.625`", "\[VeryThinSpace]", "-", 
      RowBox[{"0.625`", " ", "\[Eta]"}]}], ")"}]}]]}]], "Output", \
"PluginEmbeddedContent",
 GeneratedCell->False,
 CellAutoOverwrite->False],

Cell[BoxData[
 FractionBox[
  RowBox[{"0.625`", " ", 
   RowBox[{"(", 
    RowBox[{"1", "-", "\[Eta]"}], ")"}]}], 
  RowBox[{"5.625`", "\[VeryThinSpace]", "-", 
   RowBox[{"0.625`", " ", "\[Eta]"}]}]]], "Output", "PluginEmbeddedContent",
 GeneratedCell->False,
 CellAutoOverwrite->False],

Cell[BoxData[
 FractionBox[
  RowBox[{"0.125`", " ", 
   RowBox[{"(", 
    RowBox[{"1", "-", "\[Eta]"}], ")"}]}], 
  RowBox[{"0.875`", "\[VeryThinSpace]", "-", 
   RowBox[{"0.125`", " ", "\[Eta]"}]}]]], "Output", "PluginEmbeddedContent",
 GeneratedCell->False,
 CellAutoOverwrite->False],

Cell[BoxData[
 RowBox[{
  FractionBox[
   RowBox[{
    RowBox[{"(", 
     RowBox[{"1.75`", "\[VeryThinSpace]", "-", 
      RowBox[{"0.25`", " ", "\[Eta]"}]}], ")"}], " ", 
    RowBox[{"(", 
     RowBox[{
      RowBox[{"-", "1"}], "-", "\[Xi]"}], ")"}]}], 
   RowBox[{"4", " ", 
    RowBox[{"(", 
     RowBox[{"0.875`", "\[VeryThinSpace]", "-", 
      RowBox[{"0.125`", " ", "\[Eta]"}]}], ")"}]}]], "+", 
  FractionBox[
   RowBox[{
    RowBox[{"(", 
     RowBox[{"1", "-", "\[Eta]"}], ")"}], " ", 
    RowBox[{"(", 
     RowBox[{
      RowBox[{"-", "0.25`"}], "+", 
      RowBox[{"0.25`", " ", "\[Xi]"}]}], ")"}]}], 
   RowBox[{"4", " ", 
    RowBox[{"(", 
     RowBox[{"0.875`", "\[VeryThinSpace]", "-", 
      RowBox[{"0.125`", " ", "\[Eta]"}]}], ")"}]}]]}]], "Output", \
"PluginEmbeddedContent",
 GeneratedCell->False,
 CellAutoOverwrite->False]
}, Open  ]],

Cell[CellGroupData[{

Cell[BoxData[{
 RowBox[{"dN3dx", "=", 
  RowBox[{
   RowBox[{
    RowBox[{"D", "[", 
     RowBox[{"N3", ",", "\[Xi]"}], "]"}], "d\[Xi]dx"}], "+", 
   RowBox[{
    RowBox[{"D", "[", 
     RowBox[{"N3", ",", "\[Eta]"}], "]"}], 
    "d\[Eta]dx"}]}]}], "\[IndentingNewLine]", 
 RowBox[{"dN3dy", "=", 
  RowBox[{
   RowBox[{
    RowBox[{"D", "[", 
     RowBox[{"N3", ",", "\[Xi]"}], "]"}], "d\[Xi]dy"}], "+", 
   RowBox[{
    RowBox[{"D", "[", 
     RowBox[{"N3", ",", "\[Eta]"}], "]"}], "d\[Eta]dy"}]}]}]}], "Input", \
"PluginEmbeddedContent"],

Cell[BoxData[
 RowBox[{"0.`", "\[VeryThinSpace]", "+", 
  RowBox[{"0.1`", " ", 
   RowBox[{"(", 
    RowBox[{"1", "+", "\[Xi]"}], ")"}]}]}]], "Output", "PluginEmbeddedContent",
 GeneratedCell->False,
 CellAutoOverwrite->False],

Cell[BoxData[
 RowBox[{"0.`", "\[VeryThinSpace]", "+", 
  RowBox[{"0.1`", " ", 
   RowBox[{"(", 
    RowBox[{"1", "+", "\[Eta]"}], ")"}]}]}]], "Output", \
"PluginEmbeddedContent",
 GeneratedCell->False,
 CellAutoOverwrite->False],

Cell[BoxData[
 RowBox[{"0.`", "\[VeryThinSpace]", "+", 
  RowBox[{"0.1`", " ", 
   RowBox[{"(", 
    RowBox[{"1", "+", "\[Xi]"}], ")"}]}]}]], "Output", "PluginEmbeddedContent",
 GeneratedCell->False,
 CellAutoOverwrite->False],

Cell[BoxData[
 RowBox[{"0.`", "\[VeryThinSpace]", "+", 
  RowBox[{"0.1`", " ", 
   RowBox[{"(", 
    RowBox[{"1", "+", "\[Eta]"}], ")"}]}]}]], "Output", \
"PluginEmbeddedContent",
 GeneratedCell->False,
 CellAutoOverwrite->False],

Cell[BoxData[
 RowBox[{"0.`", "\[VeryThinSpace]", "+", 
  RowBox[{"0.1`", " ", 
   RowBox[{"(", 
    RowBox[{"1", "+", "\[Xi]"}], ")"}]}]}]], "Output", "PluginEmbeddedContent",
 GeneratedCell->False,
 CellAutoOverwrite->False],

Cell[BoxData[
 RowBox[{"0.`", "\[VeryThinSpace]", "+", 
  RowBox[{"0.1`", " ", 
   RowBox[{"(", 
    RowBox[{"1", "+", "\[Eta]"}], ")"}]}]}]], "Output", \
"PluginEmbeddedContent",
 GeneratedCell->False,
 CellAutoOverwrite->False],

Cell[BoxData[
 RowBox[{
  FractionBox[
   RowBox[{
    RowBox[{"(", 
     RowBox[{"1", "+", "\[Eta]"}], ")"}], " ", 
    RowBox[{"(", 
     RowBox[{
      RowBox[{"-", "0.25`"}], "+", 
      RowBox[{"0.25`", " ", "\[Xi]"}]}], ")"}]}], 
   RowBox[{"4", " ", 
    RowBox[{"(", 
     RowBox[{"5.625`", "\[VeryThinSpace]", "-", 
      RowBox[{"0.625`", " ", "\[Eta]"}]}], ")"}]}]], "+", 
  FractionBox[
   RowBox[{
    RowBox[{"(", 
     RowBox[{"2.25`", "\[VeryThinSpace]", "-", 
      RowBox[{"0.25`", " ", "\[Eta]"}]}], ")"}], " ", 
    RowBox[{"(", 
     RowBox[{"1", "+", "\[Xi]"}], ")"}]}], 
   RowBox[{"4", " ", 
    RowBox[{"(", 
     RowBox[{"5.625`", "\[VeryThinSpace]", "-", 
      RowBox[{"0.625`", " ", "\[Eta]"}]}], ")"}]}]]}]], "Output", \
"PluginEmbeddedContent",
 GeneratedCell->False,
 CellAutoOverwrite->False],

Cell[BoxData[
 FractionBox[
  RowBox[{"0.625`", " ", 
   RowBox[{"(", 
    RowBox[{"1", "+", "\[Eta]"}], ")"}]}], 
  RowBox[{"5.625`", "\[VeryThinSpace]", "-", 
   RowBox[{"0.625`", " ", "\[Eta]"}]}]]], "Output", "PluginEmbeddedContent",
 GeneratedCell->False,
 CellAutoOverwrite->False],

Cell[BoxData[
 FractionBox[
  RowBox[{"0.125`", " ", 
   RowBox[{"(", 
    RowBox[{"1", "+", "\[Eta]"}], ")"}]}], 
  RowBox[{"0.875`", "\[VeryThinSpace]", "-", 
   RowBox[{"0.125`", " ", "\[Eta]"}]}]]], "Output", "PluginEmbeddedContent",
 GeneratedCell->False,
 CellAutoOverwrite->False],

Cell[BoxData[
 RowBox[{
  FractionBox[
   RowBox[{
    RowBox[{"(", 
     RowBox[{"1", "+", "\[Eta]"}], ")"}], " ", 
    RowBox[{"(", 
     RowBox[{
      RowBox[{"-", "0.25`"}], "+", 
      RowBox[{"0.25`", " ", "\[Xi]"}]}], ")"}]}], 
   RowBox[{"4", " ", 
    RowBox[{"(", 
     RowBox[{"0.875`", "\[VeryThinSpace]", "-", 
      RowBox[{"0.125`", " ", "\[Eta]"}]}], ")"}]}]], "+", 
  FractionBox[
   RowBox[{
    RowBox[{"(", 
     RowBox[{"1.75`", "\[VeryThinSpace]", "-", 
      RowBox[{"0.25`", " ", "\[Eta]"}]}], ")"}], " ", 
    RowBox[{"(", 
     RowBox[{"1", "+", "\[Xi]"}], ")"}]}], 
   RowBox[{"4", " ", 
    RowBox[{"(", 
     RowBox[{"0.875`", "\[VeryThinSpace]", "-", 
      RowBox[{"0.125`", " ", "\[Eta]"}]}], ")"}]}]]}]], "Output", \
"PluginEmbeddedContent",
 GeneratedCell->False,
 CellAutoOverwrite->False]
}, Open  ]],

Cell[CellGroupData[{

Cell[BoxData[{
 RowBox[{"dN4dx", "=", 
  RowBox[{
   RowBox[{
    RowBox[{"D", "[", 
     RowBox[{"N4", ",", "\[Xi]"}], "]"}], "d\[Xi]dx"}], "+", 
   RowBox[{
    RowBox[{"D", "[", 
     RowBox[{"N4", ",", "\[Eta]"}], "]"}], 
    "d\[Eta]dx"}]}]}], "\[IndentingNewLine]", 
 RowBox[{"dN4dy", "=", 
  RowBox[{
   RowBox[{
    RowBox[{"D", "[", 
     RowBox[{"N4", ",", "\[Xi]"}], "]"}], "d\[Xi]dy"}], "+", 
   RowBox[{
    RowBox[{"D", "[", 
     RowBox[{"N4", ",", "\[Eta]"}], "]"}], "d\[Eta]dy"}]}]}]}], "Input", \
"PluginEmbeddedContent"],

Cell[BoxData[
 RowBox[{"0.`", "\[VeryThinSpace]", "+", 
  RowBox[{"0.1`", " ", 
   RowBox[{"(", 
    RowBox[{"1", "-", "\[Xi]"}], ")"}]}]}]], "Output", "PluginEmbeddedContent",
 GeneratedCell->False,
 CellAutoOverwrite->False],

Cell[BoxData[
 RowBox[{"0.`", "\[VeryThinSpace]", "+", 
  RowBox[{"0.1`", " ", 
   RowBox[{"(", 
    RowBox[{
     RowBox[{"-", "1"}], "-", "\[Eta]"}], ")"}]}]}]], "Output", \
"PluginEmbeddedContent",
 GeneratedCell->False,
 CellAutoOverwrite->False],

Cell[BoxData[
 RowBox[{"0.`", "\[VeryThinSpace]", "+", 
  RowBox[{"0.1`", " ", 
   RowBox[{"(", 
    RowBox[{"1", "-", "\[Xi]"}], ")"}]}]}]], "Output", "PluginEmbeddedContent",
 GeneratedCell->False,
 CellAutoOverwrite->False],

Cell[BoxData[
 RowBox[{"0.`", "\[VeryThinSpace]", "+", 
  RowBox[{"0.1`", " ", 
   RowBox[{"(", 
    RowBox[{
     RowBox[{"-", "1"}], "-", "\[Eta]"}], ")"}]}]}]], "Output", \
"PluginEmbeddedContent",
 GeneratedCell->False,
 CellAutoOverwrite->False],

Cell[BoxData[
 RowBox[{"0.`", "\[VeryThinSpace]", "+", 
  RowBox[{"0.1`", " ", 
   RowBox[{"(", 
    RowBox[{"1", "-", "\[Xi]"}], ")"}]}]}]], "Output", "PluginEmbeddedContent",
 GeneratedCell->False,
 CellAutoOverwrite->False],

Cell[BoxData[
 RowBox[{"0.`", "\[VeryThinSpace]", "+", 
  RowBox[{"0.1`", " ", 
   RowBox[{"(", 
    RowBox[{
     RowBox[{"-", "1"}], "-", "\[Eta]"}], ")"}]}]}]], "Output", \
"PluginEmbeddedContent",
 GeneratedCell->False,
 CellAutoOverwrite->False],

Cell[BoxData[
 RowBox[{
  FractionBox[
   RowBox[{
    RowBox[{"(", 
     RowBox[{"2.25`", "\[VeryThinSpace]", "-", 
      RowBox[{"0.25`", " ", "\[Eta]"}]}], ")"}], " ", 
    RowBox[{"(", 
     RowBox[{"1", "-", "\[Xi]"}], ")"}]}], 
   RowBox[{"4", " ", 
    RowBox[{"(", 
     RowBox[{"5.625`", "\[VeryThinSpace]", "-", 
      RowBox[{"0.625`", " ", "\[Eta]"}]}], ")"}]}]], "+", 
  FractionBox[
   RowBox[{
    RowBox[{"(", 
     RowBox[{
      RowBox[{"-", "1"}], "-", "\[Eta]"}], ")"}], " ", 
    RowBox[{"(", 
     RowBox[{
      RowBox[{"-", "0.25`"}], "+", 
      RowBox[{"0.25`", " ", "\[Xi]"}]}], ")"}]}], 
   RowBox[{"4", " ", 
    RowBox[{"(", 
     RowBox[{"5.625`", "\[VeryThinSpace]", "-", 
      RowBox[{"0.625`", " ", "\[Eta]"}]}], ")"}]}]]}]], "Output", \
"PluginEmbeddedContent",
 GeneratedCell->False,
 CellAutoOverwrite->False],

Cell[BoxData[
 FractionBox[
  RowBox[{"0.625`", " ", 
   RowBox[{"(", 
    RowBox[{
     RowBox[{"-", "1"}], "-", "\[Eta]"}], ")"}]}], 
  RowBox[{"5.625`", "\[VeryThinSpace]", "-", 
   RowBox[{"0.625`", " ", "\[Eta]"}]}]]], "Output", "PluginEmbeddedContent",
 GeneratedCell->False,
 CellAutoOverwrite->False],

Cell[BoxData[
 FractionBox[
  RowBox[{"0.125`", " ", 
   RowBox[{"(", 
    RowBox[{
     RowBox[{"-", "1"}], "-", "\[Eta]"}], ")"}]}], 
  RowBox[{"0.875`", "\[VeryThinSpace]", "-", 
   RowBox[{"0.125`", " ", "\[Eta]"}]}]]], "Output", "PluginEmbeddedContent",
 GeneratedCell->False,
 CellAutoOverwrite->False],

Cell[BoxData[
 RowBox[{
  FractionBox[
   RowBox[{
    RowBox[{"(", 
     RowBox[{"1.75`", "\[VeryThinSpace]", "-", 
      RowBox[{"0.25`", " ", "\[Eta]"}]}], ")"}], " ", 
    RowBox[{"(", 
     RowBox[{"1", "-", "\[Xi]"}], ")"}]}], 
   RowBox[{"4", " ", 
    RowBox[{"(", 
     RowBox[{"0.875`", "\[VeryThinSpace]", "-", 
      RowBox[{"0.125`", " ", "\[Eta]"}]}], ")"}]}]], "+", 
  FractionBox[
   RowBox[{
    RowBox[{"(", 
     RowBox[{
      RowBox[{"-", "1"}], "-", "\[Eta]"}], ")"}], " ", 
    RowBox[{"(", 
     RowBox[{
      RowBox[{"-", "0.25`"}], "+", 
      RowBox[{"0.25`", " ", "\[Xi]"}]}], ")"}]}], 
   RowBox[{"4", " ", 
    RowBox[{"(", 
     RowBox[{"0.875`", "\[VeryThinSpace]", "-", 
      RowBox[{"0.125`", " ", "\[Eta]"}]}], ")"}]}]]}]], "Output", \
"PluginEmbeddedContent",
 GeneratedCell->False,
 CellAutoOverwrite->False]
}, Open  ]],

Cell["Luego ensamblamos la matriz B:", "Text", "PluginEmbeddedContent"],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{
  RowBox[{"(", 
   RowBox[{"Bmat", "=", 
    RowBox[{"{", 
     RowBox[{
      RowBox[{"{", 
       RowBox[{
       "dN1dx", ",", "0", ",", "dN2dx", ",", "0", ",", "dN3dx", ",", "0", ",",
         "dN4dx", ",", "0"}], "}"}], ",", "\[IndentingNewLine]", 
      RowBox[{"{", 
       RowBox[{
       "0", ",", "dN1dy", ",", "0", ",", "dN2dy", ",", "0", ",", "dN3dy", ",",
         "0", ",", "dN4dy"}], "}"}], ",", "\[IndentingNewLine]", 
      RowBox[{"{", 
       RowBox[{
       "dN1dy", ",", "dN1dx", ",", "dN2dy", ",", "dN2dx", ",", "dN3dy", ",", 
        "dN3dx", ",", "dN4dy", ",", "dN4dx"}], "}"}]}], "}"}]}], ")"}], "//", 
  "MatrixForm"}]], "Input", "PluginEmbeddedContent"],

Cell[BoxData[
 TagBox[
  RowBox[{"(", "\[NoBreak]", GridBox[{
     {
      RowBox[{"0.`", "\[VeryThinSpace]", "+", 
       RowBox[{"0.1`", " ", 
        RowBox[{"(", 
         RowBox[{
          RowBox[{"-", "1"}], "+", "\[Eta]"}], ")"}]}]}], "0", 
      RowBox[{"0.`", "\[VeryThinSpace]", "+", 
       RowBox[{"0.1`", " ", 
        RowBox[{"(", 
         RowBox[{"1", "-", "\[Eta]"}], ")"}]}]}], "0", 
      RowBox[{"0.`", "\[VeryThinSpace]", "+", 
       RowBox[{"0.1`", " ", 
        RowBox[{"(", 
         RowBox[{"1", "+", "\[Eta]"}], ")"}]}]}], "0", 
      RowBox[{"0.`", "\[VeryThinSpace]", "+", 
       RowBox[{"0.1`", " ", 
        RowBox[{"(", 
         RowBox[{
          RowBox[{"-", "1"}], "-", "\[Eta]"}], ")"}]}]}], "0"},
     {"0", 
      RowBox[{"0.`", "\[VeryThinSpace]", "+", 
       RowBox[{"0.1`", " ", 
        RowBox[{"(", 
         RowBox[{
          RowBox[{"-", "1"}], "+", "\[Xi]"}], ")"}]}]}], "0", 
      RowBox[{"0.`", "\[VeryThinSpace]", "+", 
       RowBox[{"0.1`", " ", 
        RowBox[{"(", 
         RowBox[{
          RowBox[{"-", "1"}], "-", "\[Xi]"}], ")"}]}]}], "0", 
      RowBox[{"0.`", "\[VeryThinSpace]", "+", 
       RowBox[{"0.1`", " ", 
        RowBox[{"(", 
         RowBox[{"1", "+", "\[Xi]"}], ")"}]}]}], "0", 
      RowBox[{"0.`", "\[VeryThinSpace]", "+", 
       RowBox[{"0.1`", " ", 
        RowBox[{"(", 
         RowBox[{"1", "-", "\[Xi]"}], ")"}]}]}]},
     {
      RowBox[{"0.`", "\[VeryThinSpace]", "+", 
       RowBox[{"0.1`", " ", 
        RowBox[{"(", 
         RowBox[{
          RowBox[{"-", "1"}], "+", "\[Xi]"}], ")"}]}]}], 
      RowBox[{"0.`", "\[VeryThinSpace]", "+", 
       RowBox[{"0.1`", " ", 
        RowBox[{"(", 
         RowBox[{
          RowBox[{"-", "1"}], "+", "\[Eta]"}], ")"}]}]}], 
      RowBox[{"0.`", "\[VeryThinSpace]", "+", 
       RowBox[{"0.1`", " ", 
        RowBox[{"(", 
         RowBox[{
          RowBox[{"-", "1"}], "-", "\[Xi]"}], ")"}]}]}], 
      RowBox[{"0.`", "\[VeryThinSpace]", "+", 
       RowBox[{"0.1`", " ", 
        RowBox[{"(", 
         RowBox[{"1", "-", "\[Eta]"}], ")"}]}]}], 
      RowBox[{"0.`", "\[VeryThinSpace]", "+", 
       RowBox[{"0.1`", " ", 
        RowBox[{"(", 
         RowBox[{"1", "+", "\[Xi]"}], ")"}]}]}], 
      RowBox[{"0.`", "\[VeryThinSpace]", "+", 
       RowBox[{"0.1`", " ", 
        RowBox[{"(", 
         RowBox[{"1", "+", "\[Eta]"}], ")"}]}]}], 
      RowBox[{"0.`", "\[VeryThinSpace]", "+", 
       RowBox[{"0.1`", " ", 
        RowBox[{"(", 
         RowBox[{"1", "-", "\[Xi]"}], ")"}]}]}], 
      RowBox[{"0.`", "\[VeryThinSpace]", "+", 
       RowBox[{"0.1`", " ", 
        RowBox[{"(", 
         RowBox[{
          RowBox[{"-", "1"}], "-", "\[Eta]"}], ")"}]}]}]}
    },
    GridBoxAlignment->{
     "Columns" -> {{Center}}, "ColumnsIndexed" -> {}, "Rows" -> {{Baseline}}, 
      "RowsIndexed" -> {}},
    GridBoxSpacings->{"Columns" -> {
        Offset[0.27999999999999997`], {
         Offset[0.7]}, 
        Offset[0.27999999999999997`]}, "ColumnsIndexed" -> {}, "Rows" -> {
        Offset[0.2], {
         Offset[0.4]}, 
        Offset[0.2]}, "RowsIndexed" -> {}}], "\[NoBreak]", ")"}],
  Function[BoxForm`e$, 
   MatrixForm[BoxForm`e$]]]], "Output", "PluginEmbeddedContent",
 GeneratedCell->False,
 CellAutoOverwrite->False],

Cell[BoxData[
 TagBox[
  RowBox[{"(", "\[NoBreak]", GridBox[{
     {
      RowBox[{"0.`", "\[VeryThinSpace]", "+", 
       RowBox[{"0.1`", " ", 
        RowBox[{"(", 
         RowBox[{
          RowBox[{"-", "1"}], "+", "\[Eta]"}], ")"}]}]}], "0", 
      RowBox[{"0.`", "\[VeryThinSpace]", "+", 
       RowBox[{"0.1`", " ", 
        RowBox[{"(", 
         RowBox[{"1", "-", "\[Eta]"}], ")"}]}]}], "0", 
      RowBox[{"0.`", "\[VeryThinSpace]", "+", 
       RowBox[{"0.1`", " ", 
        RowBox[{"(", 
         RowBox[{"1", "+", "\[Eta]"}], ")"}]}]}], "0", 
      RowBox[{"0.`", "\[VeryThinSpace]", "+", 
       RowBox[{"0.1`", " ", 
        RowBox[{"(", 
         RowBox[{
          RowBox[{"-", "1"}], "-", "\[Eta]"}], ")"}]}]}], "0"},
     {"0", 
      RowBox[{"0.`", "\[VeryThinSpace]", "+", 
       RowBox[{"0.1`", " ", 
        RowBox[{"(", 
         RowBox[{
          RowBox[{"-", "1"}], "+", "\[Xi]"}], ")"}]}]}], "0", 
      RowBox[{"0.`", "\[VeryThinSpace]", "+", 
       RowBox[{"0.1`", " ", 
        RowBox[{"(", 
         RowBox[{
          RowBox[{"-", "1"}], "-", "\[Xi]"}], ")"}]}]}], "0", 
      RowBox[{"0.`", "\[VeryThinSpace]", "+", 
       RowBox[{"0.1`", " ", 
        RowBox[{"(", 
         RowBox[{"1", "+", "\[Xi]"}], ")"}]}]}], "0", 
      RowBox[{"0.`", "\[VeryThinSpace]", "+", 
       RowBox[{"0.1`", " ", 
        RowBox[{"(", 
         RowBox[{"1", "-", "\[Xi]"}], ")"}]}]}]},
     {
      RowBox[{"0.`", "\[VeryThinSpace]", "+", 
       RowBox[{"0.1`", " ", 
        RowBox[{"(", 
         RowBox[{
          RowBox[{"-", "1"}], "+", "\[Xi]"}], ")"}]}]}], 
      RowBox[{"0.`", "\[VeryThinSpace]", "+", 
       RowBox[{"0.1`", " ", 
        RowBox[{"(", 
         RowBox[{
          RowBox[{"-", "1"}], "+", "\[Eta]"}], ")"}]}]}], 
      RowBox[{"0.`", "\[VeryThinSpace]", "+", 
       RowBox[{"0.1`", " ", 
        RowBox[{"(", 
         RowBox[{
          RowBox[{"-", "1"}], "-", "\[Xi]"}], ")"}]}]}], 
      RowBox[{"0.`", "\[VeryThinSpace]", "+", 
       RowBox[{"0.1`", " ", 
        RowBox[{"(", 
         RowBox[{"1", "-", "\[Eta]"}], ")"}]}]}], 
      RowBox[{"0.`", "\[VeryThinSpace]", "+", 
       RowBox[{"0.1`", " ", 
        RowBox[{"(", 
         RowBox[{"1", "+", "\[Xi]"}], ")"}]}]}], 
      RowBox[{"0.`", "\[VeryThinSpace]", "+", 
       RowBox[{"0.1`", " ", 
        RowBox[{"(", 
         RowBox[{"1", "+", "\[Eta]"}], ")"}]}]}], 
      RowBox[{"0.`", "\[VeryThinSpace]", "+", 
       RowBox[{"0.1`", " ", 
        RowBox[{"(", 
         RowBox[{"1", "-", "\[Xi]"}], ")"}]}]}], 
      RowBox[{"0.`", "\[VeryThinSpace]", "+", 
       RowBox[{"0.1`", " ", 
        RowBox[{"(", 
         RowBox[{
          RowBox[{"-", "1"}], "-", "\[Eta]"}], ")"}]}]}]}
    },
    GridBoxAlignment->{
     "Columns" -> {{Center}}, "ColumnsIndexed" -> {}, "Rows" -> {{Baseline}}, 
      "RowsIndexed" -> {}},
    GridBoxSpacings->{"Columns" -> {
        Offset[0.27999999999999997`], {
         Offset[0.7]}, 
        Offset[0.27999999999999997`]}, "ColumnsIndexed" -> {}, "Rows" -> {
        Offset[0.2], {
         Offset[0.4]}, 
        Offset[0.2]}, "RowsIndexed" -> {}}], "\[NoBreak]", ")"}],
  Function[BoxForm`e$, 
   MatrixForm[BoxForm`e$]]]], "Output", "PluginEmbeddedContent",
 GeneratedCell->False,
 CellAutoOverwrite->False],

Cell[BoxData[
 TagBox[
  RowBox[{"(", "\[NoBreak]", GridBox[{
     {
      RowBox[{"0.`", "\[VeryThinSpace]", "+", 
       RowBox[{"0.1`", " ", 
        RowBox[{"(", 
         RowBox[{
          RowBox[{"-", "1"}], "+", "\[Eta]"}], ")"}]}]}], "0", 
      RowBox[{"0.`", "\[VeryThinSpace]", "+", 
       RowBox[{"0.1`", " ", 
        RowBox[{"(", 
         RowBox[{"1", "-", "\[Eta]"}], ")"}]}]}], "0", 
      RowBox[{"0.`", "\[VeryThinSpace]", "+", 
       RowBox[{"0.1`", " ", 
        RowBox[{"(", 
         RowBox[{"1", "+", "\[Eta]"}], ")"}]}]}], "0", 
      RowBox[{"0.`", "\[VeryThinSpace]", "+", 
       RowBox[{"0.1`", " ", 
        RowBox[{"(", 
         RowBox[{
          RowBox[{"-", "1"}], "-", "\[Eta]"}], ")"}]}]}], "0"},
     {"0", 
      RowBox[{"0.`", "\[VeryThinSpace]", "+", 
       RowBox[{"0.1`", " ", 
        RowBox[{"(", 
         RowBox[{
          RowBox[{"-", "1"}], "+", "\[Xi]"}], ")"}]}]}], "0", 
      RowBox[{"0.`", "\[VeryThinSpace]", "+", 
       RowBox[{"0.1`", " ", 
        RowBox[{"(", 
         RowBox[{
          RowBox[{"-", "1"}], "-", "\[Xi]"}], ")"}]}]}], "0", 
      RowBox[{"0.`", "\[VeryThinSpace]", "+", 
       RowBox[{"0.1`", " ", 
        RowBox[{"(", 
         RowBox[{"1", "+", "\[Xi]"}], ")"}]}]}], "0", 
      RowBox[{"0.`", "\[VeryThinSpace]", "+", 
       RowBox[{"0.1`", " ", 
        RowBox[{"(", 
         RowBox[{"1", "-", "\[Xi]"}], ")"}]}]}]},
     {
      RowBox[{"0.`", "\[VeryThinSpace]", "+", 
       RowBox[{"0.1`", " ", 
        RowBox[{"(", 
         RowBox[{
          RowBox[{"-", "1"}], "+", "\[Xi]"}], ")"}]}]}], 
      RowBox[{"0.`", "\[VeryThinSpace]", "+", 
       RowBox[{"0.1`", " ", 
        RowBox[{"(", 
         RowBox[{
          RowBox[{"-", "1"}], "+", "\[Eta]"}], ")"}]}]}], 
      RowBox[{"0.`", "\[VeryThinSpace]", "+", 
       RowBox[{"0.1`", " ", 
        RowBox[{"(", 
         RowBox[{
          RowBox[{"-", "1"}], "-", "\[Xi]"}], ")"}]}]}], 
      RowBox[{"0.`", "\[VeryThinSpace]", "+", 
       RowBox[{"0.1`", " ", 
        RowBox[{"(", 
         RowBox[{"1", "-", "\[Eta]"}], ")"}]}]}], 
      RowBox[{"0.`", "\[VeryThinSpace]", "+", 
       RowBox[{"0.1`", " ", 
        RowBox[{"(", 
         RowBox[{"1", "+", "\[Xi]"}], ")"}]}]}], 
      RowBox[{"0.`", "\[VeryThinSpace]", "+", 
       RowBox[{"0.1`", " ", 
        RowBox[{"(", 
         RowBox[{"1", "+", "\[Eta]"}], ")"}]}]}], 
      RowBox[{"0.`", "\[VeryThinSpace]", "+", 
       RowBox[{"0.1`", " ", 
        RowBox[{"(", 
         RowBox[{"1", "-", "\[Xi]"}], ")"}]}]}], 
      RowBox[{"0.`", "\[VeryThinSpace]", "+", 
       RowBox[{"0.1`", " ", 
        RowBox[{"(", 
         RowBox[{
          RowBox[{"-", "1"}], "-", "\[Eta]"}], ")"}]}]}]}
    },
    GridBoxAlignment->{
     "Columns" -> {{Center}}, "ColumnsIndexed" -> {}, "Rows" -> {{Baseline}}, 
      "RowsIndexed" -> {}},
    GridBoxSpacings->{"Columns" -> {
        Offset[0.27999999999999997`], {
         Offset[0.7]}, 
        Offset[0.27999999999999997`]}, "ColumnsIndexed" -> {}, "Rows" -> {
        Offset[0.2], {
         Offset[0.4]}, 
        Offset[0.2]}, "RowsIndexed" -> {}}], "\[NoBreak]", ")"}],
  Function[BoxForm`e$, 
   MatrixForm[BoxForm`e$]]]], "Output", "PluginEmbeddedContent",
 GeneratedCell->False,
 CellAutoOverwrite->False],

Cell[BoxData[
 TagBox[
  RowBox[{"(", "\[NoBreak]", GridBox[{
     {
      FractionBox[
       RowBox[{"0.625`", " ", 
        RowBox[{"(", 
         RowBox[{
          RowBox[{"-", "1"}], "+", "\[Eta]"}], ")"}]}], 
       RowBox[{"5.625`", "\[VeryThinSpace]", "-", 
        RowBox[{"0.625`", " ", "\[Eta]"}]}]], "0", 
      FractionBox[
       RowBox[{"0.625`", " ", 
        RowBox[{"(", 
         RowBox[{"1", "-", "\[Eta]"}], ")"}]}], 
       RowBox[{"5.625`", "\[VeryThinSpace]", "-", 
        RowBox[{"0.625`", " ", "\[Eta]"}]}]], "0", 
      FractionBox[
       RowBox[{"0.625`", " ", 
        RowBox[{"(", 
         RowBox[{"1", "+", "\[Eta]"}], ")"}]}], 
       RowBox[{"5.625`", "\[VeryThinSpace]", "-", 
        RowBox[{"0.625`", " ", "\[Eta]"}]}]], "0", 
      FractionBox[
       RowBox[{"0.625`", " ", 
        RowBox[{"(", 
         RowBox[{
          RowBox[{"-", "1"}], "-", "\[Eta]"}], ")"}]}], 
       RowBox[{"5.625`", "\[VeryThinSpace]", "-", 
        RowBox[{"0.625`", " ", "\[Eta]"}]}]], "0"},
     {"0", 
      RowBox[{
       FractionBox[
        RowBox[{
         RowBox[{"(", 
          RowBox[{
           RowBox[{"-", "1"}], "+", "\[Eta]"}], ")"}], " ", 
         RowBox[{"(", 
          RowBox[{
           RowBox[{"-", "0.25`"}], "+", 
           RowBox[{"0.25`", " ", "\[Xi]"}]}], ")"}]}], 
        RowBox[{"4", " ", 
         RowBox[{"(", 
          RowBox[{"5.625`", "\[VeryThinSpace]", "-", 
           RowBox[{"0.625`", " ", "\[Eta]"}]}], ")"}]}]], "+", 
       FractionBox[
        RowBox[{
         RowBox[{"(", 
          RowBox[{"2.25`", "\[VeryThinSpace]", "-", 
           RowBox[{"0.25`", " ", "\[Eta]"}]}], ")"}], " ", 
         RowBox[{"(", 
          RowBox[{
           RowBox[{"-", "1"}], "+", "\[Xi]"}], ")"}]}], 
        RowBox[{"4", " ", 
         RowBox[{"(", 
          RowBox[{"5.625`", "\[VeryThinSpace]", "-", 
           RowBox[{"0.625`", " ", "\[Eta]"}]}], ")"}]}]]}], "0", 
      RowBox[{
       FractionBox[
        RowBox[{
         RowBox[{"(", 
          RowBox[{"2.25`", "\[VeryThinSpace]", "-", 
           RowBox[{"0.25`", " ", "\[Eta]"}]}], ")"}], " ", 
         RowBox[{"(", 
          RowBox[{
           RowBox[{"-", "1"}], "-", "\[Xi]"}], ")"}]}], 
        RowBox[{"4", " ", 
         RowBox[{"(", 
          RowBox[{"5.625`", "\[VeryThinSpace]", "-", 
           RowBox[{"0.625`", " ", "\[Eta]"}]}], ")"}]}]], "+", 
       FractionBox[
        RowBox[{
         RowBox[{"(", 
          RowBox[{"1", "-", "\[Eta]"}], ")"}], " ", 
         RowBox[{"(", 
          RowBox[{
           RowBox[{"-", "0.25`"}], "+", 
           RowBox[{"0.25`", " ", "\[Xi]"}]}], ")"}]}], 
        RowBox[{"4", " ", 
         RowBox[{"(", 
          RowBox[{"5.625`", "\[VeryThinSpace]", "-", 
           RowBox[{"0.625`", " ", "\[Eta]"}]}], ")"}]}]]}], "0", 
      RowBox[{
       FractionBox[
        RowBox[{
         RowBox[{"(", 
          RowBox[{"1", "+", "\[Eta]"}], ")"}], " ", 
         RowBox[{"(", 
          RowBox[{
           RowBox[{"-", "0.25`"}], "+", 
           RowBox[{"0.25`", " ", "\[Xi]"}]}], ")"}]}], 
        RowBox[{"4", " ", 
         RowBox[{"(", 
          RowBox[{"5.625`", "\[VeryThinSpace]", "-", 
           RowBox[{"0.625`", " ", "\[Eta]"}]}], ")"}]}]], "+", 
       FractionBox[
        RowBox[{
         RowBox[{"(", 
          RowBox[{"2.25`", "\[VeryThinSpace]", "-", 
           RowBox[{"0.25`", " ", "\[Eta]"}]}], ")"}], " ", 
         RowBox[{"(", 
          RowBox[{"1", "+", "\[Xi]"}], ")"}]}], 
        RowBox[{"4", " ", 
         RowBox[{"(", 
          RowBox[{"5.625`", "\[VeryThinSpace]", "-", 
           RowBox[{"0.625`", " ", "\[Eta]"}]}], ")"}]}]]}], "0", 
      RowBox[{
       FractionBox[
        RowBox[{
         RowBox[{"(", 
          RowBox[{"2.25`", "\[VeryThinSpace]", "-", 
           RowBox[{"0.25`", " ", "\[Eta]"}]}], ")"}], " ", 
         RowBox[{"(", 
          RowBox[{"1", "-", "\[Xi]"}], ")"}]}], 
        RowBox[{"4", " ", 
         RowBox[{"(", 
          RowBox[{"5.625`", "\[VeryThinSpace]", "-", 
           RowBox[{"0.625`", " ", "\[Eta]"}]}], ")"}]}]], "+", 
       FractionBox[
        RowBox[{
         RowBox[{"(", 
          RowBox[{
           RowBox[{"-", "1"}], "-", "\[Eta]"}], ")"}], " ", 
         RowBox[{"(", 
          RowBox[{
           RowBox[{"-", "0.25`"}], "+", 
           RowBox[{"0.25`", " ", "\[Xi]"}]}], ")"}]}], 
        RowBox[{"4", " ", 
         RowBox[{"(", 
          RowBox[{"5.625`", "\[VeryThinSpace]", "-", 
           RowBox[{"0.625`", " ", "\[Eta]"}]}], ")"}]}]]}]},
     {
      RowBox[{
       FractionBox[
        RowBox[{
         RowBox[{"(", 
          RowBox[{
           RowBox[{"-", "1"}], "+", "\[Eta]"}], ")"}], " ", 
         RowBox[{"(", 
          RowBox[{
           RowBox[{"-", "0.25`"}], "+", 
           RowBox[{"0.25`", " ", "\[Xi]"}]}], ")"}]}], 
        RowBox[{"4", " ", 
         RowBox[{"(", 
          RowBox[{"5.625`", "\[VeryThinSpace]", "-", 
           RowBox[{"0.625`", " ", "\[Eta]"}]}], ")"}]}]], "+", 
       FractionBox[
        RowBox[{
         RowBox[{"(", 
          RowBox[{"2.25`", "\[VeryThinSpace]", "-", 
           RowBox[{"0.25`", " ", "\[Eta]"}]}], ")"}], " ", 
         RowBox[{"(", 
          RowBox[{
           RowBox[{"-", "1"}], "+", "\[Xi]"}], ")"}]}], 
        RowBox[{"4", " ", 
         RowBox[{"(", 
          RowBox[{"5.625`", "\[VeryThinSpace]", "-", 
           RowBox[{"0.625`", " ", "\[Eta]"}]}], ")"}]}]]}], 
      FractionBox[
       RowBox[{"0.625`", " ", 
        RowBox[{"(", 
         RowBox[{
          RowBox[{"-", "1"}], "+", "\[Eta]"}], ")"}]}], 
       RowBox[{"5.625`", "\[VeryThinSpace]", "-", 
        RowBox[{"0.625`", " ", "\[Eta]"}]}]], 
      RowBox[{
       FractionBox[
        RowBox[{
         RowBox[{"(", 
          RowBox[{"2.25`", "\[VeryThinSpace]", "-", 
           RowBox[{"0.25`", " ", "\[Eta]"}]}], ")"}], " ", 
         RowBox[{"(", 
          RowBox[{
           RowBox[{"-", "1"}], "-", "\[Xi]"}], ")"}]}], 
        RowBox[{"4", " ", 
         RowBox[{"(", 
          RowBox[{"5.625`", "\[VeryThinSpace]", "-", 
           RowBox[{"0.625`", " ", "\[Eta]"}]}], ")"}]}]], "+", 
       FractionBox[
        RowBox[{
         RowBox[{"(", 
          RowBox[{"1", "-", "\[Eta]"}], ")"}], " ", 
         RowBox[{"(", 
          RowBox[{
           RowBox[{"-", "0.25`"}], "+", 
           RowBox[{"0.25`", " ", "\[Xi]"}]}], ")"}]}], 
        RowBox[{"4", " ", 
         RowBox[{"(", 
          RowBox[{"5.625`", "\[VeryThinSpace]", "-", 
           RowBox[{"0.625`", " ", "\[Eta]"}]}], ")"}]}]]}], 
      FractionBox[
       RowBox[{"0.625`", " ", 
        RowBox[{"(", 
         RowBox[{"1", "-", "\[Eta]"}], ")"}]}], 
       RowBox[{"5.625`", "\[VeryThinSpace]", "-", 
        RowBox[{"0.625`", " ", "\[Eta]"}]}]], 
      RowBox[{
       FractionBox[
        RowBox[{
         RowBox[{"(", 
          RowBox[{"1", "+", "\[Eta]"}], ")"}], " ", 
         RowBox[{"(", 
          RowBox[{
           RowBox[{"-", "0.25`"}], "+", 
           RowBox[{"0.25`", " ", "\[Xi]"}]}], ")"}]}], 
        RowBox[{"4", " ", 
         RowBox[{"(", 
          RowBox[{"5.625`", "\[VeryThinSpace]", "-", 
           RowBox[{"0.625`", " ", "\[Eta]"}]}], ")"}]}]], "+", 
       FractionBox[
        RowBox[{
         RowBox[{"(", 
          RowBox[{"2.25`", "\[VeryThinSpace]", "-", 
           RowBox[{"0.25`", " ", "\[Eta]"}]}], ")"}], " ", 
         RowBox[{"(", 
          RowBox[{"1", "+", "\[Xi]"}], ")"}]}], 
        RowBox[{"4", " ", 
         RowBox[{"(", 
          RowBox[{"5.625`", "\[VeryThinSpace]", "-", 
           RowBox[{"0.625`", " ", "\[Eta]"}]}], ")"}]}]]}], 
      FractionBox[
       RowBox[{"0.625`", " ", 
        RowBox[{"(", 
         RowBox[{"1", "+", "\[Eta]"}], ")"}]}], 
       RowBox[{"5.625`", "\[VeryThinSpace]", "-", 
        RowBox[{"0.625`", " ", "\[Eta]"}]}]], 
      RowBox[{
       FractionBox[
        RowBox[{
         RowBox[{"(", 
          RowBox[{"2.25`", "\[VeryThinSpace]", "-", 
           RowBox[{"0.25`", " ", "\[Eta]"}]}], ")"}], " ", 
         RowBox[{"(", 
          RowBox[{"1", "-", "\[Xi]"}], ")"}]}], 
        RowBox[{"4", " ", 
         RowBox[{"(", 
          RowBox[{"5.625`", "\[VeryThinSpace]", "-", 
           RowBox[{"0.625`", " ", "\[Eta]"}]}], ")"}]}]], "+", 
       FractionBox[
        RowBox[{
         RowBox[{"(", 
          RowBox[{
           RowBox[{"-", "1"}], "-", "\[Eta]"}], ")"}], " ", 
         RowBox[{"(", 
          RowBox[{
           RowBox[{"-", "0.25`"}], "+", 
           RowBox[{"0.25`", " ", "\[Xi]"}]}], ")"}]}], 
        RowBox[{"4", " ", 
         RowBox[{"(", 
          RowBox[{"5.625`", "\[VeryThinSpace]", "-", 
           RowBox[{"0.625`", " ", "\[Eta]"}]}], ")"}]}]]}], 
      FractionBox[
       RowBox[{"0.625`", " ", 
        RowBox[{"(", 
         RowBox[{
          RowBox[{"-", "1"}], "-", "\[Eta]"}], ")"}]}], 
       RowBox[{"5.625`", "\[VeryThinSpace]", "-", 
        RowBox[{"0.625`", " ", "\[Eta]"}]}]]}
    },
    GridBoxAlignment->{
     "Columns" -> {{Center}}, "ColumnsIndexed" -> {}, "Rows" -> {{Baseline}}, 
      "RowsIndexed" -> {}},
    GridBoxSpacings->{"Columns" -> {
        Offset[0.27999999999999997`], {
         Offset[0.7]}, 
        Offset[0.27999999999999997`]}, "ColumnsIndexed" -> {}, "Rows" -> {
        Offset[0.2], {
         Offset[0.4]}, 
        Offset[0.2]}, "RowsIndexed" -> {}}], "\[NoBreak]", ")"}],
  Function[BoxForm`e$, 
   MatrixForm[BoxForm`e$]]]], "Output", "PluginEmbeddedContent",
 GeneratedCell->False,
 CellAutoOverwrite->False],

Cell[BoxData[
 TagBox[
  RowBox[{"(", "\[NoBreak]", GridBox[{
     {
      FractionBox[
       RowBox[{"0.125`", " ", 
        RowBox[{"(", 
         RowBox[{
          RowBox[{"-", "1"}], "+", "\[Eta]"}], ")"}]}], 
       RowBox[{"0.875`", "\[VeryThinSpace]", "-", 
        RowBox[{"0.125`", " ", "\[Eta]"}]}]], "0", 
      FractionBox[
       RowBox[{"0.125`", " ", 
        RowBox[{"(", 
         RowBox[{"1", "-", "\[Eta]"}], ")"}]}], 
       RowBox[{"0.875`", "\[VeryThinSpace]", "-", 
        RowBox[{"0.125`", " ", "\[Eta]"}]}]], "0", 
      FractionBox[
       RowBox[{"0.125`", " ", 
        RowBox[{"(", 
         RowBox[{"1", "+", "\[Eta]"}], ")"}]}], 
       RowBox[{"0.875`", "\[VeryThinSpace]", "-", 
        RowBox[{"0.125`", " ", "\[Eta]"}]}]], "0", 
      FractionBox[
       RowBox[{"0.125`", " ", 
        RowBox[{"(", 
         RowBox[{
          RowBox[{"-", "1"}], "-", "\[Eta]"}], ")"}]}], 
       RowBox[{"0.875`", "\[VeryThinSpace]", "-", 
        RowBox[{"0.125`", " ", "\[Eta]"}]}]], "0"},
     {"0", 
      RowBox[{
       FractionBox[
        RowBox[{
         RowBox[{"(", 
          RowBox[{
           RowBox[{"-", "1"}], "+", "\[Eta]"}], ")"}], " ", 
         RowBox[{"(", 
          RowBox[{
           RowBox[{"-", "0.25`"}], "+", 
           RowBox[{"0.25`", " ", "\[Xi]"}]}], ")"}]}], 
        RowBox[{"4", " ", 
         RowBox[{"(", 
          RowBox[{"0.875`", "\[VeryThinSpace]", "-", 
           RowBox[{"0.125`", " ", "\[Eta]"}]}], ")"}]}]], "+", 
       FractionBox[
        RowBox[{
         RowBox[{"(", 
          RowBox[{"1.75`", "\[VeryThinSpace]", "-", 
           RowBox[{"0.25`", " ", "\[Eta]"}]}], ")"}], " ", 
         RowBox[{"(", 
          RowBox[{
           RowBox[{"-", "1"}], "+", "\[Xi]"}], ")"}]}], 
        RowBox[{"4", " ", 
         RowBox[{"(", 
          RowBox[{"0.875`", "\[VeryThinSpace]", "-", 
           RowBox[{"0.125`", " ", "\[Eta]"}]}], ")"}]}]]}], "0", 
      RowBox[{
       FractionBox[
        RowBox[{
         RowBox[{"(", 
          RowBox[{"1.75`", "\[VeryThinSpace]", "-", 
           RowBox[{"0.25`", " ", "\[Eta]"}]}], ")"}], " ", 
         RowBox[{"(", 
          RowBox[{
           RowBox[{"-", "1"}], "-", "\[Xi]"}], ")"}]}], 
        RowBox[{"4", " ", 
         RowBox[{"(", 
          RowBox[{"0.875`", "\[VeryThinSpace]", "-", 
           RowBox[{"0.125`", " ", "\[Eta]"}]}], ")"}]}]], "+", 
       FractionBox[
        RowBox[{
         RowBox[{"(", 
          RowBox[{"1", "-", "\[Eta]"}], ")"}], " ", 
         RowBox[{"(", 
          RowBox[{
           RowBox[{"-", "0.25`"}], "+", 
           RowBox[{"0.25`", " ", "\[Xi]"}]}], ")"}]}], 
        RowBox[{"4", " ", 
         RowBox[{"(", 
          RowBox[{"0.875`", "\[VeryThinSpace]", "-", 
           RowBox[{"0.125`", " ", "\[Eta]"}]}], ")"}]}]]}], "0", 
      RowBox[{
       FractionBox[
        RowBox[{
         RowBox[{"(", 
          RowBox[{"1", "+", "\[Eta]"}], ")"}], " ", 
         RowBox[{"(", 
          RowBox[{
           RowBox[{"-", "0.25`"}], "+", 
           RowBox[{"0.25`", " ", "\[Xi]"}]}], ")"}]}], 
        RowBox[{"4", " ", 
         RowBox[{"(", 
          RowBox[{"0.875`", "\[VeryThinSpace]", "-", 
           RowBox[{"0.125`", " ", "\[Eta]"}]}], ")"}]}]], "+", 
       FractionBox[
        RowBox[{
         RowBox[{"(", 
          RowBox[{"1.75`", "\[VeryThinSpace]", "-", 
           RowBox[{"0.25`", " ", "\[Eta]"}]}], ")"}], " ", 
         RowBox[{"(", 
          RowBox[{"1", "+", "\[Xi]"}], ")"}]}], 
        RowBox[{"4", " ", 
         RowBox[{"(", 
          RowBox[{"0.875`", "\[VeryThinSpace]", "-", 
           RowBox[{"0.125`", " ", "\[Eta]"}]}], ")"}]}]]}], "0", 
      RowBox[{
       FractionBox[
        RowBox[{
         RowBox[{"(", 
          RowBox[{"1.75`", "\[VeryThinSpace]", "-", 
           RowBox[{"0.25`", " ", "\[Eta]"}]}], ")"}], " ", 
         RowBox[{"(", 
          RowBox[{"1", "-", "\[Xi]"}], ")"}]}], 
        RowBox[{"4", " ", 
         RowBox[{"(", 
          RowBox[{"0.875`", "\[VeryThinSpace]", "-", 
           RowBox[{"0.125`", " ", "\[Eta]"}]}], ")"}]}]], "+", 
       FractionBox[
        RowBox[{
         RowBox[{"(", 
          RowBox[{
           RowBox[{"-", "1"}], "-", "\[Eta]"}], ")"}], " ", 
         RowBox[{"(", 
          RowBox[{
           RowBox[{"-", "0.25`"}], "+", 
           RowBox[{"0.25`", " ", "\[Xi]"}]}], ")"}]}], 
        RowBox[{"4", " ", 
         RowBox[{"(", 
          RowBox[{"0.875`", "\[VeryThinSpace]", "-", 
           RowBox[{"0.125`", " ", "\[Eta]"}]}], ")"}]}]]}]},
     {
      RowBox[{
       FractionBox[
        RowBox[{
         RowBox[{"(", 
          RowBox[{
           RowBox[{"-", "1"}], "+", "\[Eta]"}], ")"}], " ", 
         RowBox[{"(", 
          RowBox[{
           RowBox[{"-", "0.25`"}], "+", 
           RowBox[{"0.25`", " ", "\[Xi]"}]}], ")"}]}], 
        RowBox[{"4", " ", 
         RowBox[{"(", 
          RowBox[{"0.875`", "\[VeryThinSpace]", "-", 
           RowBox[{"0.125`", " ", "\[Eta]"}]}], ")"}]}]], "+", 
       FractionBox[
        RowBox[{
         RowBox[{"(", 
          RowBox[{"1.75`", "\[VeryThinSpace]", "-", 
           RowBox[{"0.25`", " ", "\[Eta]"}]}], ")"}], " ", 
         RowBox[{"(", 
          RowBox[{
           RowBox[{"-", "1"}], "+", "\[Xi]"}], ")"}]}], 
        RowBox[{"4", " ", 
         RowBox[{"(", 
          RowBox[{"0.875`", "\[VeryThinSpace]", "-", 
           RowBox[{"0.125`", " ", "\[Eta]"}]}], ")"}]}]]}], 
      FractionBox[
       RowBox[{"0.125`", " ", 
        RowBox[{"(", 
         RowBox[{
          RowBox[{"-", "1"}], "+", "\[Eta]"}], ")"}]}], 
       RowBox[{"0.875`", "\[VeryThinSpace]", "-", 
        RowBox[{"0.125`", " ", "\[Eta]"}]}]], 
      RowBox[{
       FractionBox[
        RowBox[{
         RowBox[{"(", 
          RowBox[{"1.75`", "\[VeryThinSpace]", "-", 
           RowBox[{"0.25`", " ", "\[Eta]"}]}], ")"}], " ", 
         RowBox[{"(", 
          RowBox[{
           RowBox[{"-", "1"}], "-", "\[Xi]"}], ")"}]}], 
        RowBox[{"4", " ", 
         RowBox[{"(", 
          RowBox[{"0.875`", "\[VeryThinSpace]", "-", 
           RowBox[{"0.125`", " ", "\[Eta]"}]}], ")"}]}]], "+", 
       FractionBox[
        RowBox[{
         RowBox[{"(", 
          RowBox[{"1", "-", "\[Eta]"}], ")"}], " ", 
         RowBox[{"(", 
          RowBox[{
           RowBox[{"-", "0.25`"}], "+", 
           RowBox[{"0.25`", " ", "\[Xi]"}]}], ")"}]}], 
        RowBox[{"4", " ", 
         RowBox[{"(", 
          RowBox[{"0.875`", "\[VeryThinSpace]", "-", 
           RowBox[{"0.125`", " ", "\[Eta]"}]}], ")"}]}]]}], 
      FractionBox[
       RowBox[{"0.125`", " ", 
        RowBox[{"(", 
         RowBox[{"1", "-", "\[Eta]"}], ")"}]}], 
       RowBox[{"0.875`", "\[VeryThinSpace]", "-", 
        RowBox[{"0.125`", " ", "\[Eta]"}]}]], 
      RowBox[{
       FractionBox[
        RowBox[{
         RowBox[{"(", 
          RowBox[{"1", "+", "\[Eta]"}], ")"}], " ", 
         RowBox[{"(", 
          RowBox[{
           RowBox[{"-", "0.25`"}], "+", 
           RowBox[{"0.25`", " ", "\[Xi]"}]}], ")"}]}], 
        RowBox[{"4", " ", 
         RowBox[{"(", 
          RowBox[{"0.875`", "\[VeryThinSpace]", "-", 
           RowBox[{"0.125`", " ", "\[Eta]"}]}], ")"}]}]], "+", 
       FractionBox[
        RowBox[{
         RowBox[{"(", 
          RowBox[{"1.75`", "\[VeryThinSpace]", "-", 
           RowBox[{"0.25`", " ", "\[Eta]"}]}], ")"}], " ", 
         RowBox[{"(", 
          RowBox[{"1", "+", "\[Xi]"}], ")"}]}], 
        RowBox[{"4", " ", 
         RowBox[{"(", 
          RowBox[{"0.875`", "\[VeryThinSpace]", "-", 
           RowBox[{"0.125`", " ", "\[Eta]"}]}], ")"}]}]]}], 
      FractionBox[
       RowBox[{"0.125`", " ", 
        RowBox[{"(", 
         RowBox[{"1", "+", "\[Eta]"}], ")"}]}], 
       RowBox[{"0.875`", "\[VeryThinSpace]", "-", 
        RowBox[{"0.125`", " ", "\[Eta]"}]}]], 
      RowBox[{
       FractionBox[
        RowBox[{
         RowBox[{"(", 
          RowBox[{"1.75`", "\[VeryThinSpace]", "-", 
           RowBox[{"0.25`", " ", "\[Eta]"}]}], ")"}], " ", 
         RowBox[{"(", 
          RowBox[{"1", "-", "\[Xi]"}], ")"}]}], 
        RowBox[{"4", " ", 
         RowBox[{"(", 
          RowBox[{"0.875`", "\[VeryThinSpace]", "-", 
           RowBox[{"0.125`", " ", "\[Eta]"}]}], ")"}]}]], "+", 
       FractionBox[
        RowBox[{
         RowBox[{"(", 
          RowBox[{
           RowBox[{"-", "1"}], "-", "\[Eta]"}], ")"}], " ", 
         RowBox[{"(", 
          RowBox[{
           RowBox[{"-", "0.25`"}], "+", 
           RowBox[{"0.25`", " ", "\[Xi]"}]}], ")"}]}], 
        RowBox[{"4", " ", 
         RowBox[{"(", 
          RowBox[{"0.875`", "\[VeryThinSpace]", "-", 
           RowBox[{"0.125`", " ", "\[Eta]"}]}], ")"}]}]]}], 
      FractionBox[
       RowBox[{"0.125`", " ", 
        RowBox[{"(", 
         RowBox[{
          RowBox[{"-", "1"}], "-", "\[Eta]"}], ")"}]}], 
       RowBox[{"0.875`", "\[VeryThinSpace]", "-", 
        RowBox[{"0.125`", " ", "\[Eta]"}]}]]}
    },
    GridBoxAlignment->{
     "Columns" -> {{Center}}, "ColumnsIndexed" -> {}, "Rows" -> {{Baseline}}, 
      "RowsIndexed" -> {}},
    GridBoxSpacings->{"Columns" -> {
        Offset[0.27999999999999997`], {
         Offset[0.7]}, 
        Offset[0.27999999999999997`]}, "ColumnsIndexed" -> {}, "Rows" -> {
        Offset[0.2], {
         Offset[0.4]}, 
        Offset[0.2]}, "RowsIndexed" -> {}}], "\[NoBreak]", ")"}],
  Function[BoxForm`e$, 
   MatrixForm[BoxForm`e$]]]], "Output", "PluginEmbeddedContent",
 GeneratedCell->False,
 CellAutoOverwrite->False]
}, Open  ]],

Cell[CellGroupData[{

Cell[BoxData[{
 RowBox[{
  RowBox[{"\[Nu]", "=", "0.33"}], ";"}], "\[IndentingNewLine]", 
 RowBox[{"Y", "=", 
  RowBox[{"10.0", "*", 
   RowBox[{"10", "^", "6"}]}]}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{"(", 
   RowBox[{"Cmat", "=", 
    RowBox[{
     RowBox[{"Y", "/", 
      RowBox[{"(", 
       RowBox[{"1", "-", 
        RowBox[{"\[Nu]", "^", "2"}]}], ")"}]}], 
     RowBox[{"{", 
      RowBox[{
       RowBox[{"{", 
        RowBox[{"1", ",", "\[Nu]", ",", "0"}], "}"}], ",", 
       RowBox[{"{", 
        RowBox[{"\[Nu]", ",", "1", ",", "0"}], "}"}], ",", 
       RowBox[{"{", 
        RowBox[{"0", ",", "0", ",", 
         RowBox[{
          RowBox[{"(", 
           RowBox[{"1", "-", "\[Nu]"}], ")"}], "/", "2"}]}], "}"}]}], 
      "}"}]}]}], ")"}], "//", "MatrixForm"}]}], "Input", \
"PluginEmbeddedContent"],

Cell[BoxData[
 TagBox[
  RowBox[{"(", "\[NoBreak]", GridBox[{
     {"1.1222085063404782`*^7", "3.703288070923578`*^6", "0.`"},
     {"3.703288070923578`*^6", "1.1222085063404782`*^7", "0.`"},
     {"0.`", "0.`", "3.7593984962406014`*^6"}
    },
    GridBoxAlignment->{
     "Columns" -> {{Center}}, "ColumnsIndexed" -> {}, "Rows" -> {{Baseline}}, 
      "RowsIndexed" -> {}},
    GridBoxSpacings->{"Columns" -> {
        Offset[0.27999999999999997`], {
         Offset[0.7]}, 
        Offset[0.27999999999999997`]}, "ColumnsIndexed" -> {}, "Rows" -> {
        Offset[0.2], {
         Offset[0.4]}, 
        Offset[0.2]}, "RowsIndexed" -> {}}], "\[NoBreak]", ")"}],
  Function[BoxForm`e$, 
   MatrixForm[BoxForm`e$]]]], "Output", "PluginEmbeddedContent",
 GeneratedCell->False,
 CellAutoOverwrite->False],

Cell[BoxData["1.`*^7"], "Output", "PluginEmbeddedContent",
 GeneratedCell->False,
 CellAutoOverwrite->False],

Cell[BoxData[
 TagBox[
  RowBox[{"(", "\[NoBreak]", GridBox[{
     {"1.1222085063404782`*^7", "3.703288070923578`*^6", "0.`"},
     {"3.703288070923578`*^6", "1.1222085063404782`*^7", "0.`"},
     {"0.`", "0.`", "3.7593984962406014`*^6"}
    },
    GridBoxAlignment->{
     "Columns" -> {{Center}}, "ColumnsIndexed" -> {}, "Rows" -> {{Baseline}}, 
      "RowsIndexed" -> {}},
    GridBoxSpacings->{"Columns" -> {
        Offset[0.27999999999999997`], {
         Offset[0.7]}, 
        Offset[0.27999999999999997`]}, "ColumnsIndexed" -> {}, "Rows" -> {
        Offset[0.2], {
         Offset[0.4]}, 
        Offset[0.2]}, "RowsIndexed" -> {}}], "\[NoBreak]", ")"}],
  Function[BoxForm`e$, 
   MatrixForm[BoxForm`e$]]]], "Output", "PluginEmbeddedContent",
 GeneratedCell->False,
 CellAutoOverwrite->False],

Cell[BoxData["1.`*^7"], "Output", "PluginEmbeddedContent",
 GeneratedCell->False,
 CellAutoOverwrite->False],

Cell[BoxData[
 TagBox[
  RowBox[{"(", "\[NoBreak]", GridBox[{
     {"1.1222085063404782`*^7", "3.703288070923578`*^6", "0.`"},
     {"3.703288070923578`*^6", "1.1222085063404782`*^7", "0.`"},
     {"0.`", "0.`", "3.7593984962406014`*^6"}
    },
    GridBoxAlignment->{
     "Columns" -> {{Center}}, "ColumnsIndexed" -> {}, "Rows" -> {{Baseline}}, 
      "RowsIndexed" -> {}},
    GridBoxSpacings->{"Columns" -> {
        Offset[0.27999999999999997`], {
         Offset[0.7]}, 
        Offset[0.27999999999999997`]}, "ColumnsIndexed" -> {}, "Rows" -> {
        Offset[0.2], {
         Offset[0.4]}, 
        Offset[0.2]}, "RowsIndexed" -> {}}], "\[NoBreak]", ")"}],
  Function[BoxForm`e$, 
   MatrixForm[BoxForm`e$]]]], "Output", "PluginEmbeddedContent",
 GeneratedCell->False,
 CellAutoOverwrite->False],

Cell[BoxData["1.`*^7"], "Output", "PluginEmbeddedContent",
 GeneratedCell->False,
 CellAutoOverwrite->False],

Cell[BoxData[
 TagBox[
  RowBox[{"(", "\[NoBreak]", GridBox[{
     {"1.1222085063404782`*^7", "3.703288070923578`*^6", "0.`"},
     {"3.703288070923578`*^6", "1.1222085063404782`*^7", "0.`"},
     {"0.`", "0.`", "3.7593984962406014`*^6"}
    },
    GridBoxAlignment->{
     "Columns" -> {{Center}}, "ColumnsIndexed" -> {}, "Rows" -> {{Baseline}}, 
      "RowsIndexed" -> {}},
    GridBoxSpacings->{"Columns" -> {
        Offset[0.27999999999999997`], {
         Offset[0.7]}, 
        Offset[0.27999999999999997`]}, "ColumnsIndexed" -> {}, "Rows" -> {
        Offset[0.2], {
         Offset[0.4]}, 
        Offset[0.2]}, "RowsIndexed" -> {}}], "\[NoBreak]", ")"}],
  Function[BoxForm`e$, 
   MatrixForm[BoxForm`e$]]]], "Output", "PluginEmbeddedContent",
 GeneratedCell->False,
 CellAutoOverwrite->False],

Cell[BoxData["1.`*^7"], "Output", "PluginEmbeddedContent",
 GeneratedCell->False,
 CellAutoOverwrite->False],

Cell[BoxData["1.`*^7"], "Output", "PluginEmbeddedContent",
 GeneratedCell->False,
 CellAutoOverwrite->False],

Cell[BoxData[
 TagBox[
  RowBox[{"(", "\[NoBreak]", GridBox[{
     {"1.1222085063404782`*^7", "3.703288070923578`*^6", "0.`"},
     {"3.703288070923578`*^6", "1.1222085063404782`*^7", "0.`"},
     {"0.`", "0.`", "3.7593984962406014`*^6"}
    },
    GridBoxAlignment->{
     "Columns" -> {{Center}}, "ColumnsIndexed" -> {}, "Rows" -> {{Baseline}}, 
      "RowsIndexed" -> {}},
    GridBoxSpacings->{"Columns" -> {
        Offset[0.27999999999999997`], {
         Offset[0.7]}, 
        Offset[0.27999999999999997`]}, "ColumnsIndexed" -> {}, "Rows" -> {
        Offset[0.2], {
         Offset[0.4]}, 
        Offset[0.2]}, "RowsIndexed" -> {}}], "\[NoBreak]", ")"}],
  Function[BoxForm`e$, 
   MatrixForm[BoxForm`e$]]]], "Output", "PluginEmbeddedContent",
 GeneratedCell->False,
 CellAutoOverwrite->False]
}, Open  ]],

Cell["\<\
Ahora utilizamos cuadratura gaussinana para calcular la integral\
\>", "Text", "PluginEmbeddedContent"],

Cell[BoxData[{
 RowBox[{
  RowBox[{"q1", "=", 
   RowBox[{"-", "0.5773502691"}]}], ";", 
  RowBox[{"q2", "=", 
   RowBox[{"-", "q1"}]}], ";", 
  RowBox[{"w1", "=", "1.0"}], ";", 
  RowBox[{"w2", "=", "1.0"}], ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{"intg", "=", 
   RowBox[{
    RowBox[{"(", 
     RowBox[{
      RowBox[{"Transpose", "[", "Bmat", "]"}], ".", "Cmat", ".", "Bmat"}], 
     ")"}], "jac"}]}], ";"}]}], "Input", "PluginEmbeddedContent"],

Cell["Ahora la matriz de rigidez es:", "Text", "PluginEmbeddedContent"],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{
  RowBox[{"(", 
   RowBox[{"ke", "=", 
    RowBox[{
     RowBox[{"(", 
      RowBox[{
       RowBox[{"intg", " ", "w1", " ", "w1"}], "/.", 
       RowBox[{"{", 
        RowBox[{
         RowBox[{"\[Xi]", "\[Rule]", "q1"}], ",", 
         RowBox[{"\[Eta]", "\[Rule]", "q1"}]}], "}"}]}], ")"}], "+", 
     RowBox[{"(", 
      RowBox[{
       RowBox[{"intg", " ", "w1", " ", "w2"}], "/.", 
       RowBox[{"{", 
        RowBox[{
         RowBox[{"\[Xi]", "\[Rule]", "q1"}], ",", 
         RowBox[{"\[Eta]", "\[Rule]", "q2"}]}], "}"}]}], ")"}], "+", 
     "\[IndentingNewLine]", 
     RowBox[{"(", 
      RowBox[{
       RowBox[{"intg", " ", "w2", " ", "w1"}], "/.", 
       RowBox[{"{", 
        RowBox[{
         RowBox[{"\[Xi]", "\[Rule]", "q2"}], ",", 
         RowBox[{"\[Eta]", "\[Rule]", "q1"}]}], "}"}]}], ")"}], "+", 
     RowBox[{"(", 
      RowBox[{
       RowBox[{"intg", " ", "w2", " ", "w2"}], "/.", 
       RowBox[{"{", 
        RowBox[{
         RowBox[{"\[Xi]", "\[Rule]", "q2"}], ",", 
         RowBox[{"\[Eta]", "\[Rule]", "q2"}]}], "}"}]}], ")"}]}]}], ")"}], "//",
   "MatrixForm"}]], "Input", "PluginEmbeddedContent"],

Cell[BoxData[
 TagBox[
  RowBox[{"(", "\[NoBreak]", GridBox[{
     {"4.993827852827517`*^6", "1.8656716417910454`*^6", 
      RowBox[{"-", "3.114128604707216`*^6"}], 
      RowBox[{"-", "14027.606329255854`"}], 
      RowBox[{"-", "2.496913926995177`*^6"}], 
      RowBox[{"-", "1.8656716417910454`*^6"}], "617214.6788748753`", 
      "14027.606329255868`"},
     {"1.8656716417910454`*^6", "4.993827852827517`*^6", 
      "14027.606329255854`", "617214.6788748753`", 
      RowBox[{"-", "1.8656716417910454`*^6"}], 
      RowBox[{"-", "2.496913926995177`*^6"}], 
      RowBox[{"-", "14027.606329255868`"}], 
      RowBox[{"-", "3.114128604707216`*^6"}]},
     {
      RowBox[{"-", "3.114128604707216`*^6"}], "14027.606329255854`", 
      "4.993827852827517`*^6", 
      RowBox[{"-", "1.8656716417910454`*^6"}], "617214.6788748753`", 
      RowBox[{"-", "14027.606329255868`"}], 
      RowBox[{"-", "2.496913926995177`*^6"}], "1.8656716417910454`*^6"},
     {
      RowBox[{"-", "14027.606329255854`"}], "617214.6788748753`", 
      RowBox[{"-", "1.8656716417910454`*^6"}], "4.993827852827517`*^6", 
      "14027.606329255868`", 
      RowBox[{"-", "3.114128604707216`*^6"}], "1.8656716417910454`*^6", 
      RowBox[{"-", "2.496913926995177`*^6"}]},
     {
      RowBox[{"-", "2.496913926995177`*^6"}], 
      RowBox[{"-", "1.8656716417910454`*^6"}], "617214.6788748753`", 
      "14027.606329255868`", "4.993827852827517`*^6", 
      "1.8656716417910454`*^6", 
      RowBox[{"-", "3.114128604707216`*^6"}], 
      RowBox[{"-", "14027.606329255854`"}]},
     {
      RowBox[{"-", "1.8656716417910454`*^6"}], 
      RowBox[{"-", "2.496913926995177`*^6"}], 
      RowBox[{"-", "14027.606329255868`"}], 
      RowBox[{"-", "3.114128604707216`*^6"}], "1.8656716417910454`*^6", 
      "4.993827852827517`*^6", "14027.606329255854`", "617214.6788748753`"},
     {"617214.6788748753`", 
      RowBox[{"-", "14027.606329255868`"}], 
      RowBox[{"-", "2.496913926995177`*^6"}], "1.8656716417910454`*^6", 
      RowBox[{"-", "3.114128604707216`*^6"}], "14027.606329255854`", 
      "4.993827852827517`*^6", 
      RowBox[{"-", "1.8656716417910454`*^6"}]},
     {"14027.606329255868`", 
      RowBox[{"-", "3.114128604707216`*^6"}], "1.8656716417910454`*^6", 
      RowBox[{"-", "2.496913926995177`*^6"}], 
      RowBox[{"-", "14027.606329255854`"}], "617214.6788748753`", 
      RowBox[{"-", "1.8656716417910454`*^6"}], "4.993827852827517`*^6"}
    },
    GridBoxAlignment->{
     "Columns" -> {{Center}}, "ColumnsIndexed" -> {}, "Rows" -> {{Baseline}}, 
      "RowsIndexed" -> {}},
    GridBoxSpacings->{"Columns" -> {
        Offset[0.27999999999999997`], {
         Offset[0.7]}, 
        Offset[0.27999999999999997`]}, "ColumnsIndexed" -> {}, "Rows" -> {
        Offset[0.2], {
         Offset[0.4]}, 
        Offset[0.2]}, "RowsIndexed" -> {}}], "\[NoBreak]", ")"}],
  Function[BoxForm`e$, 
   MatrixForm[BoxForm`e$]]]], "Output", "PluginEmbeddedContent",
 GeneratedCell->False,
 CellAutoOverwrite->False],

Cell[BoxData[
 TagBox[
  RowBox[{"(", "\[NoBreak]", GridBox[{
     {"4.993827852827517`*^6", "1.8656716417910454`*^6", 
      RowBox[{"-", "3.114128604707216`*^6"}], 
      RowBox[{"-", "14027.606329255854`"}], 
      RowBox[{"-", "2.496913926995177`*^6"}], 
      RowBox[{"-", "1.8656716417910454`*^6"}], "617214.6788748753`", 
      "14027.606329255868`"},
     {"1.8656716417910454`*^6", "4.993827852827517`*^6", 
      "14027.606329255854`", "617214.6788748753`", 
      RowBox[{"-", "1.8656716417910454`*^6"}], 
      RowBox[{"-", "2.496913926995177`*^6"}], 
      RowBox[{"-", "14027.606329255868`"}], 
      RowBox[{"-", "3.114128604707216`*^6"}]},
     {
      RowBox[{"-", "3.114128604707216`*^6"}], "14027.606329255854`", 
      "4.993827852827517`*^6", 
      RowBox[{"-", "1.8656716417910454`*^6"}], "617214.6788748753`", 
      RowBox[{"-", "14027.606329255868`"}], 
      RowBox[{"-", "2.496913926995177`*^6"}], "1.8656716417910454`*^6"},
     {
      RowBox[{"-", "14027.606329255854`"}], "617214.6788748753`", 
      RowBox[{"-", "1.8656716417910454`*^6"}], "4.993827852827517`*^6", 
      "14027.606329255868`", 
      RowBox[{"-", "3.114128604707216`*^6"}], "1.8656716417910454`*^6", 
      RowBox[{"-", "2.496913926995177`*^6"}]},
     {
      RowBox[{"-", "2.496913926995177`*^6"}], 
      RowBox[{"-", "1.8656716417910454`*^6"}], "617214.6788748753`", 
      "14027.606329255868`", "4.993827852827517`*^6", 
      "1.8656716417910454`*^6", 
      RowBox[{"-", "3.114128604707216`*^6"}], 
      RowBox[{"-", "14027.606329255854`"}]},
     {
      RowBox[{"-", "1.8656716417910454`*^6"}], 
      RowBox[{"-", "2.496913926995177`*^6"}], 
      RowBox[{"-", "14027.606329255868`"}], 
      RowBox[{"-", "3.114128604707216`*^6"}], "1.8656716417910454`*^6", 
      "4.993827852827517`*^6", "14027.606329255854`", "617214.6788748753`"},
     {"617214.6788748753`", 
      RowBox[{"-", "14027.606329255868`"}], 
      RowBox[{"-", "2.496913926995177`*^6"}], "1.8656716417910454`*^6", 
      RowBox[{"-", "3.114128604707216`*^6"}], "14027.606329255854`", 
      "4.993827852827517`*^6", 
      RowBox[{"-", "1.8656716417910454`*^6"}]},
     {"14027.606329255868`", 
      RowBox[{"-", "3.114128604707216`*^6"}], "1.8656716417910454`*^6", 
      RowBox[{"-", "2.496913926995177`*^6"}], 
      RowBox[{"-", "14027.606329255854`"}], "617214.6788748753`", 
      RowBox[{"-", "1.8656716417910454`*^6"}], "4.993827852827517`*^6"}
    },
    GridBoxAlignment->{
     "Columns" -> {{Center}}, "ColumnsIndexed" -> {}, "Rows" -> {{Baseline}}, 
      "RowsIndexed" -> {}},
    GridBoxSpacings->{"Columns" -> {
        Offset[0.27999999999999997`], {
         Offset[0.7]}, 
        Offset[0.27999999999999997`]}, "ColumnsIndexed" -> {}, "Rows" -> {
        Offset[0.2], {
         Offset[0.4]}, 
        Offset[0.2]}, "RowsIndexed" -> {}}], "\[NoBreak]", ")"}],
  Function[BoxForm`e$, 
   MatrixForm[BoxForm`e$]]]], "Output", "PluginEmbeddedContent",
 GeneratedCell->False,
 CellAutoOverwrite->False],

Cell[BoxData[
 TagBox[
  RowBox[{"(", "\[NoBreak]", GridBox[{
     {"4.993827852827517`*^6", "1.8656716417910454`*^6", 
      RowBox[{"-", "3.114128604707216`*^6"}], 
      RowBox[{"-", "14027.606329255854`"}], 
      RowBox[{"-", "2.496913926995177`*^6"}], 
      RowBox[{"-", "1.8656716417910454`*^6"}], "617214.6788748753`", 
      "14027.606329255868`"},
     {"1.8656716417910454`*^6", "4.993827852827517`*^6", 
      "14027.606329255854`", "617214.6788748753`", 
      RowBox[{"-", "1.8656716417910454`*^6"}], 
      RowBox[{"-", "2.496913926995177`*^6"}], 
      RowBox[{"-", "14027.606329255868`"}], 
      RowBox[{"-", "3.114128604707216`*^6"}]},
     {
      RowBox[{"-", "3.114128604707216`*^6"}], "14027.606329255854`", 
      "4.993827852827517`*^6", 
      RowBox[{"-", "1.8656716417910454`*^6"}], "617214.6788748753`", 
      RowBox[{"-", "14027.606329255868`"}], 
      RowBox[{"-", "2.496913926995177`*^6"}], "1.8656716417910454`*^6"},
     {
      RowBox[{"-", "14027.606329255854`"}], "617214.6788748753`", 
      RowBox[{"-", "1.8656716417910454`*^6"}], "4.993827852827517`*^6", 
      "14027.606329255868`", 
      RowBox[{"-", "3.114128604707216`*^6"}], "1.8656716417910454`*^6", 
      RowBox[{"-", "2.496913926995177`*^6"}]},
     {
      RowBox[{"-", "2.496913926995177`*^6"}], 
      RowBox[{"-", "1.8656716417910454`*^6"}], "617214.6788748753`", 
      "14027.606329255868`", "4.993827852827517`*^6", 
      "1.8656716417910454`*^6", 
      RowBox[{"-", "3.114128604707216`*^6"}], 
      RowBox[{"-", "14027.606329255854`"}]},
     {
      RowBox[{"-", "1.8656716417910454`*^6"}], 
      RowBox[{"-", "2.496913926995177`*^6"}], 
      RowBox[{"-", "14027.606329255868`"}], 
      RowBox[{"-", "3.114128604707216`*^6"}], "1.8656716417910454`*^6", 
      "4.993827852827517`*^6", "14027.606329255854`", "617214.6788748753`"},
     {"617214.6788748753`", 
      RowBox[{"-", "14027.606329255868`"}], 
      RowBox[{"-", "2.496913926995177`*^6"}], "1.8656716417910454`*^6", 
      RowBox[{"-", "3.114128604707216`*^6"}], "14027.606329255854`", 
      "4.993827852827517`*^6", 
      RowBox[{"-", "1.8656716417910454`*^6"}]},
     {"14027.606329255868`", 
      RowBox[{"-", "3.114128604707216`*^6"}], "1.8656716417910454`*^6", 
      RowBox[{"-", "2.496913926995177`*^6"}], 
      RowBox[{"-", "14027.606329255854`"}], "617214.6788748753`", 
      RowBox[{"-", "1.8656716417910454`*^6"}], "4.993827852827517`*^6"}
    },
    GridBoxAlignment->{
     "Columns" -> {{Center}}, "ColumnsIndexed" -> {}, "Rows" -> {{Baseline}}, 
      "RowsIndexed" -> {}},
    GridBoxSpacings->{"Columns" -> {
        Offset[0.27999999999999997`], {
         Offset[0.7]}, 
        Offset[0.27999999999999997`]}, "ColumnsIndexed" -> {}, "Rows" -> {
        Offset[0.2], {
         Offset[0.4]}, 
        Offset[0.2]}, "RowsIndexed" -> {}}], "\[NoBreak]", ")"}],
  Function[BoxForm`e$, 
   MatrixForm[BoxForm`e$]]]], "Output", "PluginEmbeddedContent",
 GeneratedCell->False,
 CellAutoOverwrite->False]
}, Open  ]],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{
  RowBox[{"(", 
   RowBox[{"keb", "=", 
    RowBox[{"{", 
     RowBox[{
      RowBox[{"ke", "[", 
       RowBox[{"[", 
        RowBox[{"1", ",", 
         RowBox[{"{", 
          RowBox[{"3", ",", "4", ",", "5", ",", "6"}], "}"}]}], "]"}], "]"}], 
      ",", 
      RowBox[{"ke", "[", 
       RowBox[{"[", 
        RowBox[{"2", ",", 
         RowBox[{"{", 
          RowBox[{"3", ",", "4", ",", "5", ",", "6"}], "}"}]}], "]"}], "]"}], 
      ",", 
      RowBox[{"ke", "[", 
       RowBox[{"[", 
        RowBox[{"3", ",", 
         RowBox[{"{", 
          RowBox[{"3", ",", "4", ",", "5", ",", "6"}], "}"}]}], "]"}], "]"}], 
      ",", "\[IndentingNewLine]", 
      RowBox[{"ke", "[", 
       RowBox[{"[", 
        RowBox[{"4", ",", 
         RowBox[{"{", 
          RowBox[{"3", ",", "4", ",", "5", ",", "6"}], "}"}]}], "]"}], "]"}], 
      ",", 
      RowBox[{"ke", "[", 
       RowBox[{"[", 
        RowBox[{"5", ",", 
         RowBox[{"{", 
          RowBox[{"3", ",", "4", ",", "5", ",", "6"}], "}"}]}], "]"}], "]"}], 
      ",", 
      RowBox[{"ke", "[", 
       RowBox[{"[", 
        RowBox[{"6", ",", 
         RowBox[{"{", 
          RowBox[{"3", ",", "4", ",", "5", ",", "6"}], "}"}]}], "]"}], "]"}], 
      ",", "\[IndentingNewLine]", 
      RowBox[{"ke", "[", 
       RowBox[{"[", 
        RowBox[{"7", ",", 
         RowBox[{"{", 
          RowBox[{"3", ",", "4", ",", "5", ",", "6"}], "}"}]}], "]"}], "]"}], 
      ",", 
      RowBox[{"ke", "[", 
       RowBox[{"[", 
        RowBox[{"8", ",", 
         RowBox[{"{", 
          RowBox[{"3", ",", "4", ",", "5", ",", "6"}], "}"}]}], "]"}], 
       "]"}]}], "}"}]}], ")"}], "//", "MatrixForm"}]], "Input", \
"PluginEmbeddedContent"],

Cell[BoxData[
 TagBox[
  RowBox[{"(", "\[NoBreak]", GridBox[{
     {
      RowBox[{"-", "3.114128604707216`*^6"}], 
      RowBox[{"-", "14027.606329255854`"}], 
      RowBox[{"-", "2.496913926995177`*^6"}], 
      RowBox[{"-", "1.8656716417910454`*^6"}]},
     {"14027.606329255854`", "617214.6788748753`", 
      RowBox[{"-", "1.8656716417910454`*^6"}], 
      RowBox[{"-", "2.496913926995177`*^6"}]},
     {"4.993827852827517`*^6", 
      RowBox[{"-", "1.8656716417910454`*^6"}], "617214.6788748753`", 
      RowBox[{"-", "14027.606329255868`"}]},
     {
      RowBox[{"-", "1.8656716417910454`*^6"}], "4.993827852827517`*^6", 
      "14027.606329255868`", 
      RowBox[{"-", "3.114128604707216`*^6"}]},
     {"617214.6788748753`", "14027.606329255868`", "4.993827852827517`*^6", 
      "1.8656716417910454`*^6"},
     {
      RowBox[{"-", "14027.606329255868`"}], 
      RowBox[{"-", "3.114128604707216`*^6"}], "1.8656716417910454`*^6", 
      "4.993827852827517`*^6"},
     {
      RowBox[{"-", "2.496913926995177`*^6"}], "1.8656716417910454`*^6", 
      RowBox[{"-", "3.114128604707216`*^6"}], "14027.606329255854`"},
     {"1.8656716417910454`*^6", 
      RowBox[{"-", "2.496913926995177`*^6"}], 
      RowBox[{"-", "14027.606329255854`"}], "617214.6788748753`"}
    },
    GridBoxAlignment->{
     "Columns" -> {{Center}}, "ColumnsIndexed" -> {}, "Rows" -> {{Baseline}}, 
      "RowsIndexed" -> {}},
    GridBoxSpacings->{"Columns" -> {
        Offset[0.27999999999999997`], {
         Offset[0.7]}, 
        Offset[0.27999999999999997`]}, "ColumnsIndexed" -> {}, "Rows" -> {
        Offset[0.2], {
         Offset[0.4]}, 
        Offset[0.2]}, "RowsIndexed" -> {}}], "\[NoBreak]", ")"}],
  Function[BoxForm`e$, 
   MatrixForm[BoxForm`e$]]]], "Output", "PluginEmbeddedContent",
 GeneratedCell->False,
 CellAutoOverwrite->False],

Cell[BoxData[
 TagBox[
  RowBox[{"(", "\[NoBreak]", GridBox[{
     {
      RowBox[{"-", "3.114128604707216`*^6"}], 
      RowBox[{"-", "14027.606329255854`"}], 
      RowBox[{"-", "2.496913926995177`*^6"}], 
      RowBox[{"-", "1.8656716417910454`*^6"}]},
     {"14027.606329255854`", "617214.6788748753`", 
      RowBox[{"-", "1.8656716417910454`*^6"}], 
      RowBox[{"-", "2.496913926995177`*^6"}]},
     {"4.993827852827517`*^6", 
      RowBox[{"-", "1.8656716417910454`*^6"}], "617214.6788748753`", 
      RowBox[{"-", "14027.606329255868`"}]},
     {
      RowBox[{"-", "1.8656716417910454`*^6"}], "4.993827852827517`*^6", 
      "14027.606329255868`", 
      RowBox[{"-", "3.114128604707216`*^6"}]},
     {"617214.6788748753`", "14027.606329255868`", "4.993827852827517`*^6", 
      "1.8656716417910454`*^6"},
     {
      RowBox[{"-", "14027.606329255868`"}], 
      RowBox[{"-", "3.114128604707216`*^6"}], "1.8656716417910454`*^6", 
      "4.993827852827517`*^6"},
     {
      RowBox[{"-", "2.496913926995177`*^6"}], "1.8656716417910454`*^6", 
      RowBox[{"-", "3.114128604707216`*^6"}], "14027.606329255854`"},
     {"1.8656716417910454`*^6", 
      RowBox[{"-", "2.496913926995177`*^6"}], 
      RowBox[{"-", "14027.606329255854`"}], "617214.6788748753`"}
    },
    GridBoxAlignment->{
     "Columns" -> {{Center}}, "ColumnsIndexed" -> {}, "Rows" -> {{Baseline}}, 
      "RowsIndexed" -> {}},
    GridBoxSpacings->{"Columns" -> {
        Offset[0.27999999999999997`], {
         Offset[0.7]}, 
        Offset[0.27999999999999997`]}, "ColumnsIndexed" -> {}, "Rows" -> {
        Offset[0.2], {
         Offset[0.4]}, 
        Offset[0.2]}, "RowsIndexed" -> {}}], "\[NoBreak]", ")"}],
  Function[BoxForm`e$, 
   MatrixForm[BoxForm`e$]]]], "Output", "PluginEmbeddedContent",
 GeneratedCell->False,
 CellAutoOverwrite->False],

Cell[BoxData[
 TagBox[
  RowBox[{"(", "\[NoBreak]", GridBox[{
     {
      RowBox[{"-", "3.114128604707216`*^6"}], 
      RowBox[{"-", "14027.606329255854`"}], 
      RowBox[{"-", "2.496913926995177`*^6"}], 
      RowBox[{"-", "1.8656716417910454`*^6"}]},
     {"14027.606329255854`", "617214.6788748753`", 
      RowBox[{"-", "1.8656716417910454`*^6"}], 
      RowBox[{"-", "2.496913926995177`*^6"}]},
     {"4.993827852827517`*^6", 
      RowBox[{"-", "1.8656716417910454`*^6"}], "617214.6788748753`", 
      RowBox[{"-", "14027.606329255868`"}]},
     {
      RowBox[{"-", "1.8656716417910454`*^6"}], "4.993827852827517`*^6", 
      "14027.606329255868`", 
      RowBox[{"-", "3.114128604707216`*^6"}]},
     {"617214.6788748753`", "14027.606329255868`", "4.993827852827517`*^6", 
      "1.8656716417910454`*^6"},
     {
      RowBox[{"-", "14027.606329255868`"}], 
      RowBox[{"-", "3.114128604707216`*^6"}], "1.8656716417910454`*^6", 
      "4.993827852827517`*^6"},
     {
      RowBox[{"-", "2.496913926995177`*^6"}], "1.8656716417910454`*^6", 
      RowBox[{"-", "3.114128604707216`*^6"}], "14027.606329255854`"},
     {"1.8656716417910454`*^6", 
      RowBox[{"-", "2.496913926995177`*^6"}], 
      RowBox[{"-", "14027.606329255854`"}], "617214.6788748753`"}
    },
    GridBoxAlignment->{
     "Columns" -> {{Center}}, "ColumnsIndexed" -> {}, "Rows" -> {{Baseline}}, 
      "RowsIndexed" -> {}},
    GridBoxSpacings->{"Columns" -> {
        Offset[0.27999999999999997`], {
         Offset[0.7]}, 
        Offset[0.27999999999999997`]}, "ColumnsIndexed" -> {}, "Rows" -> {
        Offset[0.2], {
         Offset[0.4]}, 
        Offset[0.2]}, "RowsIndexed" -> {}}], "\[NoBreak]", ")"}],
  Function[BoxForm`e$, 
   MatrixForm[BoxForm`e$]]]], "Output", "PluginEmbeddedContent",
 GeneratedCell->False,
 CellAutoOverwrite->False],

Cell[BoxData[
 TagBox[
  RowBox[{"(", "\[NoBreak]", GridBox[{
     {
      RowBox[{"-", "3.332681029989641`*^6"}], "248091.54993834515`", 
      RowBox[{"-", "2.4720120275168777`*^6"}], 
      RowBox[{"-", "1.7269026766777849`*^6"}]},
     {"276146.76259685704`", "497342.40651592216`", 
      RowBox[{"-", "1.7269026766777853`*^6"}], 
      RowBox[{"-", "1.8490935619548007`*^6"}]},
     {"5.212380278109942`*^6", 
      RowBox[{"-", "2.1277907980586467`*^6"}], "968252.629020637`", 
      RowBox[{"-", "152796.57144251576`"}]},
     {
      RowBox[{"-", "2.1277907980586467`*^6"}], "5.1137001251864685`*^6", 
      RowBox[{"-", "124741.35878400384`"}], 
      RowBox[{"-", "2.639740463407112`*^6"}]},
     {"968252.629020637`", 
      RowBox[{"-", "124741.3587840039`"}], "5.427547528728134`*^6", 
      "1.5727149377348586`*^6"},
     {
      RowBox[{"-", "152796.57144251565`"}], 
      RowBox[{"-", "2.639740463407112`*^6"}], "1.5727149377348586`*^6", 
      "4.527091133068788`*^6"},
     {
      RowBox[{"-", "2.8479518771409383`*^6"}], "2.004440606904305`*^6", 
      RowBox[{"-", "3.923788130231893`*^6"}], "306984.3103854423`"},
     {"2.004440606904305`*^6", 
      RowBox[{"-", "2.971302068295279`*^6"}], "278929.09772693063`", 
      RowBox[{"-", "38257.10770687522`"}]}
    },
    GridBoxAlignment->{
     "Columns" -> {{Center}}, "ColumnsIndexed" -> {}, "Rows" -> {{Baseline}}, 
      "RowsIndexed" -> {}},
    GridBoxSpacings->{"Columns" -> {
        Offset[0.27999999999999997`], {
         Offset[0.7]}, 
        Offset[0.27999999999999997`]}, "ColumnsIndexed" -> {}, "Rows" -> {
        Offset[0.2], {
         Offset[0.4]}, 
        Offset[0.2]}, "RowsIndexed" -> {}}], "\[NoBreak]", ")"}],
  Function[BoxForm`e$, 
   MatrixForm[BoxForm`e$]]]], "Output", "PluginEmbeddedContent",
 GeneratedCell->False,
 CellAutoOverwrite->False],

Cell[BoxData[
 TagBox[
  RowBox[{"(", "\[NoBreak]", GridBox[{
     {"1.3954585915952616`*^6", "318215.28876367444`", 
      RowBox[{"-", "1.8512597179075119`*^6"}], 
      RowBox[{"-", "1.686771621317937`*^6"}]},
     {"346270.50142218603`", "6.813573496995181`*^6", 
      RowBox[{"-", "1.686771621317937`*^6"}], 
      RowBox[{"-", "4.1002885469979504`*^6"}]},
     {"6.123338400885941`*^6", 
      RowBox[{"-", "2.197914536883975`*^6"}], 
      RowBox[{"-", "3.787838026453391`*^6"}], 
      RowBox[{"-", "192927.62680236387`"}]},
     {
      RowBox[{"-", "2.1979145368839754`*^6"}], "1.5630596629814386`*^7", 
      RowBox[{"-", "164872.41414385207`"}], 
      RowBox[{"-", "1.2732839048109222`*^7"}]},
     {
      RowBox[{"-", "3.787838026453391`*^6"}], 
      RowBox[{"-", "164872.414143852`"}], "5.041098964385018`*^6", 
      "1.4823144551581736`*^6"},
     {
      RowBox[{"-", "192927.62680236375`"}], 
      RowBox[{"-", "1.2732839048109222`*^7"}], "1.482314455158174`*^6", 
      "1.1992642615150008`*^7"},
     {
      RowBox[{"-", "3.7309589660278126`*^6"}], "2.044571662264153`*^6", 
      "597998.7799758854`", "397384.7929621269`"},
     {"2.044571662264153`*^6", 
      RowBox[{"-", "9.711331078700343`*^6"}], "369329.580303615`", 
      "4.840484979957167`*^6"}
    },
    GridBoxAlignment->{
     "Columns" -> {{Center}}, "ColumnsIndexed" -> {}, "Rows" -> {{Baseline}}, 
      "RowsIndexed" -> {}},
    GridBoxSpacings->{"Columns" -> {
        Offset[0.27999999999999997`], {
         Offset[0.7]}, 
        Offset[0.27999999999999997`]}, "ColumnsIndexed" -> {}, "Rows" -> {
        Offset[0.2], {
         Offset[0.4]}, 
        Offset[0.2]}, "RowsIndexed" -> {}}], "\[NoBreak]", ")"}],
  Function[BoxForm`e$, 
   MatrixForm[BoxForm`e$]]]], "Output", "PluginEmbeddedContent",
 GeneratedCell->False,
 CellAutoOverwrite->False]
}, Open  ]],

Cell[CellGroupData[{

Cell[BoxData[{
 RowBox[{
  RowBox[{"ue", "=", 
   RowBox[{"{", 
    RowBox[{
    "0", ",", "u2", ",", "0.1", ",", "0", ",", "0.1", ",", "0", ",", "0", ",",
      "u8"}], "}"}]}], ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{"fe", "=", 
   RowBox[{"ke", ".", "ue"}]}], ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{"eq1", "=", 
   RowBox[{
    RowBox[{"fe", "[", 
     RowBox[{"[", "1", "]"}], "]"}], "==", "f1"}]}], ";", 
  RowBox[{"eq2", "=", 
   RowBox[{
    RowBox[{"fe", "[", 
     RowBox[{"[", "2", "]"}], "]"}], "==", "0"}]}], ";", 
  RowBox[{"eq3", "=", 
   RowBox[{
    RowBox[{"fe", "[", 
     RowBox[{"[", "3", "]"}], "]"}], "==", "f3"}]}], ";", 
  RowBox[{"eq4", "=", 
   RowBox[{
    RowBox[{"fe", "[", 
     RowBox[{"[", "4", "]"}], "]"}], "==", "f4"}]}], ";", 
  RowBox[{"eq5", "=", 
   RowBox[{
    RowBox[{"fe", "[", 
     RowBox[{"[", "5", "]"}], "]"}], "==", "f5"}]}], ";", 
  RowBox[{"eq6", "=", 
   RowBox[{
    RowBox[{"fe", "[", 
     RowBox[{"[", "6", "]"}], "]"}], "==", "f6"}]}], ";", 
  RowBox[{"eq7", "=", 
   RowBox[{
    RowBox[{"fe", "[", 
     RowBox[{"[", "7", "]"}], "]"}], "==", "f7"}]}], ";", 
  RowBox[{"eq8", "=", 
   RowBox[{
    RowBox[{"fe", "[", 
     RowBox[{"[", "8", "]"}], "]"}], "==", "0"}]}], 
  ";"}], "\[IndentingNewLine]", 
 RowBox[{"Solve", "[", 
  RowBox[{
   RowBox[{
   "eq1", "&&", "eq2", "&&", "eq3", "&&", "eq4", "&&", "eq5", "&&", "eq6", "&&",
     " ", "eq7", " ", "&&", " ", "eq8"}], ",", 
   RowBox[{"{", 
    RowBox[{
    "f1", ",", "f3", ",", "f4", ",", "f5", ",", "f6", ",", "f7", ",", "u2", 
     ",", "u8"}], "}"}]}], "]"}]}], "Input", "PluginEmbeddedContent"],

Cell[BoxData[
 RowBox[{"{", 
  RowBox[{"{", 
   RowBox[{
    RowBox[{"f1", "\[Rule]", 
     RowBox[{"-", "518817.5727715671`"}]}], ",", 
    RowBox[{"f3", "\[Rule]", "518817.572771567`"}], ",", 
    RowBox[{"f4", "\[Rule]", 
     RowBox[{"-", "114045.89558525337`"}]}], ",", 
    RowBox[{"f5", "\[Rule]", "518817.5727715671`"}], ",", 
    RowBox[{"f6", "\[Rule]", "114045.89558525337`"}], ",", 
    RowBox[{"f7", "\[Rule]", 
     RowBox[{"-", "518817.572771567`"}]}], ",", 
    RowBox[{"u2", "\[Rule]", "0.022837370244398074`"}], ",", 
    RowBox[{"u8", "\[Rule]", 
     RowBox[{"-", "0.022837370244398084`"}]}]}], "}"}], "}"}]], "Output", \
"PluginEmbeddedContent"]
}, Open  ]]
}, Open  ]]
}, Open  ]],

Cell[CellGroupData[{

Cell["\<\
Implementaci\[OAcute]n de un elemento Q4:\
\>", "Section", "PluginEmbeddedContent"],

Cell["Definiendo las coordenadas de los nodos de la malla:", "Text", "PluginEmbeddedContent"],

Cell[BoxData[{
 RowBox[{"ClearAll", "[", "\"\<Global`*\>\"", "]"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{
  "dir", "=", 
   "\"\</Users/jorgedelacruz/Dropbox/Google Drive/Mathematica\>\""}], 
  ";"}], "\[IndentingNewLine]", 
 RowBox[{"Get", "[", 
  RowBox[{"\"\<EQ4.m\>\"", ",", 
   RowBox[{"Path", "\[Rule]", "dir"}]}], "]"}]}], "Input", \
"PluginEmbeddedContent"],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{"MallaQ4", "[", 
  RowBox[{"4", ",", "3", ",", "10.", ",", "5."}], "]"}]], "Input", \
"PluginEmbeddedContent",
 CellGroupingRules->{GroupTogetherGrouping, 10000.}],

Cell[BoxData[
 GraphicsBox[
  GraphicsComplexBox[{{0., 0.}, {3.333333333333333, 0.}, {6.666666666666666, 
   0.}, {10., 0.}, {0., 2.5}, {3.333333333333333, 2.5}, {6.666666666666666, 
   2.5}, {10., 2.5}, {0., 5.}, {3.333333333333333, 5.}, {6.666666666666666, 
   5.}, {10., 
   5.}}, {LineBox[{{1, 2, 6, 5, 1}, {2, 3, 7, 6, 2}, {3, 4, 8, 7, 3}, {5, 6, 
     10, 9, 5}, {6, 7, 11, 10, 6}, {7, 8, 12, 11, 7}}], 
    {RGBColor[1, 0, 0], PointSize[Large], 
     PointBox[{{0., 0.}, {3.333333333333333, 0.}, {6.666666666666666, 0.}, {
      10., 0.}, {0., 2.5}, {3.333333333333333, 2.5}, {6.666666666666666, 
      2.5}, {10., 2.5}, {0., 5.}, {3.333333333333333, 5.}, {6.666666666666666,
       5.}, {10., 5.}}]}}]]], "Output", "PluginEmbeddedContent",
 CellGroupingRules->{GroupTogetherGrouping, 10000.}],

Cell["\<\
C\[AAcute]lculo de las matrices de rigidez elementales:\
\>", "Text", "PluginEmbeddedContent",
 CellGroupingRules->{GroupTogetherGrouping, 10000.}],

Cell[BoxData[{
 RowBox[{
  RowBox[{"p0", "=", 
   RowBox[{"{", 
    RowBox[{
     RowBox[{"{", 
      RowBox[{"0.", ",", "0."}], "}"}], ",", 
     RowBox[{"{", 
      RowBox[{"3.", ",", "1."}], "}"}], ",", 
     RowBox[{"{", 
      RowBox[{"5", ",", "0.5"}], "}"}], ",", 
     RowBox[{"{", 
      RowBox[{"10.", ",", "0."}], "}"}], ",", 
     RowBox[{"{", 
      RowBox[{"0.", ",", "5."}], "}"}], ",", 
     RowBox[{"{", 
      RowBox[{"3.", ",", "4."}], "}"}], ",", 
     RowBox[{"{", 
      RowBox[{"6", ",", "4.5"}], "}"}], ",", 
     RowBox[{"{", 
      RowBox[{"10.", ",", "5."}], "}"}], ",", 
     RowBox[{"{", 
      RowBox[{"0.0", ",", "9.0"}], "}"}], ",", 
     RowBox[{"{", 
      RowBox[{"4.0", ",", "8.5"}], "}"}], ",", 
     RowBox[{"{", 
      RowBox[{"6.0", ",", "9.0"}], "}"}], ",", 
     RowBox[{"{", 
      RowBox[{"10.0", ",", "9.0"}], "}"}]}], "}"}]}], 
  ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{"GlobalMK", "[", "p0", "]"}], ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{"GlobalK", "//", "MatrixForm"}], ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{"GlobalK", "-", 
   RowBox[{"Transpose", "[", "GlobalK", "]"}]}], ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{
   RowBox[{"MBDBe", "[", 
    RowBox[{"[", "1", "]"}], "]"}], "//", "MatrixForm"}], ";"}]}], "Input", \
"PluginEmbeddedContent",
 CellGroupingRules->{GroupTogetherGrouping, 10000.}],

Cell[BoxData[{
 RowBox[{
  RowBox[{"GlobalKAmp", "=", 
   RowBox[{"Array", "[", 
    RowBox[{
     RowBox[{"0", "&"}], ",", 
     RowBox[{"{", 
      RowBox[{
       RowBox[{"2", " ", "NNodos"}], ",", 
       RowBox[{"4", " ", "NNodos"}]}], "}"}]}], "]"}]}], 
  ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{
   RowBox[{"GlobalKAmp", "[", 
    RowBox[{"[", 
     RowBox[{"{", 
      RowBox[{"1", ",", "12"}], "}"}], "]"}], "]"}], "//", "MatrixForm"}], 
  ";"}], "\[IndentingNewLine]", 
 RowBox[{"For", "[", 
  RowBox[{
   RowBox[{"i", "=", "1"}], ",", 
   RowBox[{"i", "<", 
    RowBox[{"(", 
     RowBox[{
      RowBox[{"2", " ", "NNodos"}], "+", "1"}], ")"}]}], ",", 
   RowBox[{"i", "++"}], ",", 
   RowBox[{"(", "\[IndentingNewLine]", 
    RowBox[{"For", "[", 
     RowBox[{
      RowBox[{"j", "=", "1"}], ",", 
      RowBox[{"j", "<", 
       RowBox[{"(", 
        RowBox[{
         RowBox[{"2", " ", "NNodos"}], "+", "1"}], ")"}]}], ",", 
      RowBox[{"j", "++"}], ",", 
      RowBox[{"(", "\[IndentingNewLine]", 
       RowBox[{
        RowBox[{"GlobalKAmp", "[", 
         RowBox[{"[", 
          RowBox[{"i", ",", "j"}], "]"}], "]"}], "=", 
        RowBox[{"GlobalK", "[", 
         RowBox[{"[", 
          RowBox[{"i", ",", "j"}], "]"}], "]"}]}], "\[IndentingNewLine]", 
       ")"}]}], "]"}], "\[IndentingNewLine]", ")"}]}], 
  "]"}], "\[IndentingNewLine]", 
 RowBox[{"For", "[", 
  RowBox[{
   RowBox[{"i", "=", "1"}], ",", 
   RowBox[{"i", "<", 
    RowBox[{"(", 
     RowBox[{
      RowBox[{"2", " ", "NNodos"}], "+", "1"}], ")"}]}], ",", 
   RowBox[{"i", "++"}], ",", 
   RowBox[{"(", "\[IndentingNewLine]", 
    RowBox[{
     RowBox[{
      RowBox[{"GlobalKAmp", "[", 
       RowBox[{"[", 
        RowBox[{"i", ",", 
         RowBox[{"i", "+", 
          RowBox[{"2", " ", "NNodos"}]}]}], "]"}], "]"}], "=", "1"}], ";"}], 
    "\[IndentingNewLine]", ")"}]}], "]"}], "\n", 
 RowBox[{
  RowBox[{"GlobalKAmp", "//", "MatrixForm"}], ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{"Dimensions", "[", "GlobalKAmp", "]"}], 
  ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{"VecIn", "=", 
   RowBox[{"Array", "[", 
    RowBox[{
     RowBox[{"0", "&"}], ",", 
     RowBox[{"{", 
      RowBox[{"4", " ", "NNodos"}], "}"}]}], "]"}]}], 
  ";"}], "\[IndentingNewLine]", 
 RowBox[{"DeltaU", "=", 
  RowBox[{"Array", "[", 
   RowBox[{
    RowBox[{"0", "&"}], ",", 
    RowBox[{"{", " ", 
     RowBox[{"NNodos", ",", "2"}], "}"}]}], "]"}]}]}], "Input", \
"PluginEmbeddedContent",
 CellGroupingRules->{GroupTogetherGrouping, 10000.}],

Cell[BoxData[
 RowBox[{"{", 
  RowBox[{
   RowBox[{"{", 
    RowBox[{"0", ",", "0"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"0", ",", "0"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"0", ",", "0"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"0", ",", "0"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"0", ",", "0"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"0", ",", "0"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"0", ",", "0"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"0", ",", "0"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"0", ",", "0"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"0", ",", "0"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"0", ",", "0"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"0", ",", "0"}], "}"}]}], "}"}]], "Output", \
"PluginEmbeddedContent",
 CellGroupingRules->{GroupTogetherGrouping, 10000.}],

Cell["\<\
Ahora se especifican las condiones de contorno del problema dentro del vector \
VecIn:\
\>", "Text", "PluginEmbeddedContent",
 CellGroupingRules->{GroupTogetherGrouping, 10000.}],

Cell[BoxData[{
 RowBox[{
  RowBox[{
   RowBox[{"VecIn", "[", 
    RowBox[{"[", "7", "]"}], "]"}], "=", "4.0"}], ";", 
  RowBox[{
   RowBox[{"VecIn", "[", 
    RowBox[{"[", "15", "]"}], "]"}], "=", "4.0"}], ";", 
  RowBox[{
   RowBox[{"VecIn", "[", 
    RowBox[{"[", "23", "]"}], "]"}], "=", "4.0"}], 
  ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{
   RowBox[{"VecIn", "[", 
    RowBox[{"[", "3", "]"}], "]"}], "=", "u2x"}], ";", 
  RowBox[{
   RowBox[{"VecIn", "[", 
    RowBox[{"[", "4", "]"}], "]"}], "=", "u2y"}], ";", 
  RowBox[{
   RowBox[{"VecIn", "[", 
    RowBox[{"[", "5", "]"}], "]"}], "=", "u3x"}], ";", 
  RowBox[{
   RowBox[{"VecIn", "[", 
    RowBox[{"[", "6", "]"}], "]"}], "=", "u3y"}], ";", 
  RowBox[{
   RowBox[{"VecIn", "[", 
    RowBox[{"[", "11", "]"}], "]"}], "=", "u6x"}], ";", 
  RowBox[{
   RowBox[{"VecIn", "[", 
    RowBox[{"[", "12", "]"}], "]"}], "=", "u6y"}], ";", 
  RowBox[{
   RowBox[{"VecIn", "[", 
    RowBox[{"[", "13", "]"}], "]"}], "=", "u7x"}], ";", 
  RowBox[{
   RowBox[{"VecIn", "[", 
    RowBox[{"[", "14", "]"}], "]"}], "=", "u7y"}], ";", 
  RowBox[{
   RowBox[{"VecIn", "[", 
    RowBox[{"[", "19", "]"}], "]"}], "=", "u10x"}], ";", 
  RowBox[{
   RowBox[{"VecIn", "[", 
    RowBox[{"[", "20", "]"}], "]"}], "=", "u10y"}], ";", 
  RowBox[{
   RowBox[{"VecIn", "[", 
    RowBox[{"[", "21", "]"}], "]"}], "=", "u11x"}], ";", 
  RowBox[{
   RowBox[{"VecIn", "[", 
    RowBox[{"[", "22", "]"}], "]"}], "=", "u11y"}], ";"}]}], "Input", \
"PluginEmbeddedContent",
 CellGroupingRules->{GroupTogetherGrouping, 10000.}],

Cell["Ahora se especifican las fuerzas:", "Text", "PluginEmbeddedContent",
 CellGroupingRules->{GroupTogetherGrouping, 10000.}],

Cell[BoxData[
 RowBox[{
  RowBox[{
   RowBox[{"VecIn", "[", 
    RowBox[{"[", "25", "]"}], "]"}], "=", "f1x"}], ";", 
  RowBox[{
   RowBox[{"VecIn", "[", 
    RowBox[{"[", "26", "]"}], "]"}], "=", "f1y"}], ";", 
  RowBox[{
   RowBox[{"VecIn", "[", 
    RowBox[{"[", "31", "]"}], "]"}], "=", "f4x"}], ";", 
  RowBox[{
   RowBox[{"VecIn", "[", 
    RowBox[{"[", "32", "]"}], "]"}], "=", "f4y"}], ";", 
  RowBox[{
   RowBox[{"VecIn", "[", 
    RowBox[{"[", "33", "]"}], "]"}], "=", "f5x"}], ";", 
  RowBox[{
   RowBox[{"VecIn", "[", 
    RowBox[{"[", "34", "]"}], "]"}], "=", "f5y"}], ";", 
  RowBox[{
   RowBox[{"VecIn", "[", 
    RowBox[{"[", "39", "]"}], "]"}], "=", "f8x"}], ";", 
  RowBox[{
   RowBox[{"VecIn", "[", 
    RowBox[{"[", "40", "]"}], "]"}], "=", "f8y"}], ";", 
  RowBox[{
   RowBox[{"VecIn", "[", 
    RowBox[{"[", "41", "]"}], "]"}], "=", "f9x"}], ";", 
  RowBox[{
   RowBox[{"VecIn", "[", 
    RowBox[{"[", "42", "]"}], "]"}], "=", "f9y"}], ";", 
  RowBox[{
   RowBox[{"VecIn", "[", 
    RowBox[{"[", "47", "]"}], "]"}], "=", "f12x"}], ";", 
  RowBox[{
   RowBox[{"VecIn", "[", 
    RowBox[{"[", "48", "]"}], "]"}], "=", "f12y"}], ";"}]], "Input", \
"PluginEmbeddedContent",
 CellGroupingRules->{GroupTogetherGrouping, 10000.}],

Cell[BoxData[
 RowBox[{"\[IndentingNewLine]", 
  RowBox[{
   RowBox[{
    RowBox[{"Solucion", "=", 
     RowBox[{"GlobalKAmp", ".", "VecIn"}]}], ";"}], "\[IndentingNewLine]", 
   RowBox[{
    RowBox[{"Dimensions", "[", "Solucion", "]"}], ";"}], 
   "\[IndentingNewLine]", 
   RowBox[{"Solucion", "=", 
    RowBox[{"NSolve", "[", 
     RowBox[{
      RowBox[{"{", 
       RowBox[{
        RowBox[{
         RowBox[{"Solucion", "[", 
          RowBox[{"[", "2", "]"}], "]"}], "\[Equal]", "0"}], " ", ",", " ", 
        RowBox[{
         RowBox[{"Solucion", "[", 
          RowBox[{"[", "1", "]"}], "]"}], "\[Equal]", "0"}], " ", ",", 
        RowBox[{
         RowBox[{"Solucion", "[", 
          RowBox[{"[", "3", "]"}], "]"}], "\[Equal]", "0"}], ",", 
        RowBox[{
         RowBox[{"Solucion", "[", 
          RowBox[{"[", "4", "]"}], "]"}], "\[Equal]", "0"}], ",", 
        RowBox[{
         RowBox[{"Solucion", "[", 
          RowBox[{"[", "5", "]"}], "]"}], "\[Equal]", "0"}], 
        "\[IndentingNewLine]", ",", 
        RowBox[{
         RowBox[{"Solucion", "[", 
          RowBox[{"[", "6", "]"}], "]"}], "\[Equal]", "0"}], ",", 
        RowBox[{
         RowBox[{"Solucion", "[", 
          RowBox[{"[", "7", "]"}], "]"}], "\[Equal]", "0"}], " ", ",", 
        RowBox[{
         RowBox[{"Solucion", "[", 
          RowBox[{"[", "9", "]"}], "]"}], "\[Equal]", "0"}], " ", ",", 
        RowBox[{
         RowBox[{"Solucion", "[", 
          RowBox[{"[", "8", "]"}], "]"}], "\[Equal]", "0"}], " ", ",", 
        RowBox[{
         RowBox[{"Solucion", "[", 
          RowBox[{"[", "10", "]"}], "]"}], "\[Equal]", "0"}], 
        "\[IndentingNewLine]", ",", 
        RowBox[{
         RowBox[{"Solucion", "[", 
          RowBox[{"[", "11", "]"}], "]"}], "\[Equal]", "0"}], ",", 
        RowBox[{
         RowBox[{"Solucion", "[", 
          RowBox[{"[", "12", "]"}], "]"}], "\[Equal]", "0"}], ",", 
        RowBox[{
         RowBox[{"Solucion", "[", 
          RowBox[{"[", "13", "]"}], "]"}], "\[Equal]", "0"}], ",", 
        RowBox[{
         RowBox[{"Solucion", "[", 
          RowBox[{"[", "14", "]"}], "]"}], "\[Equal]", "0"}], " ", ",", 
        RowBox[{
         RowBox[{"Solucion", "[", 
          RowBox[{"[", "15", "]"}], "]"}], "\[Equal]", "0"}], " ", ",", 
        RowBox[{
         RowBox[{"Solucion", "[", 
          RowBox[{"[", "16", "]"}], "]"}], "\[Equal]", "0"}], "  ", ",", 
        RowBox[{
         RowBox[{"Solucion", "[", 
          RowBox[{"[", "17", "]"}], "]"}], "\[Equal]", "0"}], " ", ",", 
        RowBox[{
         RowBox[{"Solucion", "[", 
          RowBox[{"[", "18", "]"}], "]"}], "\[Equal]", "0"}], " ", ",", 
        RowBox[{
         RowBox[{"Solucion", "[", 
          RowBox[{"[", "19", "]"}], "]"}], "\[Equal]", "0"}], ",", 
        RowBox[{
         RowBox[{"Solucion", "[", 
          RowBox[{"[", "20", "]"}], "]"}], "\[Equal]", "0"}], ",", 
        RowBox[{
         RowBox[{"Solucion", "[", 
          RowBox[{"[", "21", "]"}], "]"}], "\[Equal]", "0"}], ",", 
        RowBox[{
         RowBox[{"Solucion", "[", 
          RowBox[{"[", "22", "]"}], "]"}], "\[Equal]", "0"}], " ", ",", 
        RowBox[{
         RowBox[{"Solucion", "[", 
          RowBox[{"[", "23", "]"}], "]"}], "\[Equal]", "0"}], " ", ",", 
        RowBox[{
         RowBox[{"Solucion", "[", 
          RowBox[{"[", "24", "]"}], "]"}], "\[Equal]", "0"}]}], " ", " ", " ",
        " ", "}"}], ",", 
      RowBox[{"{", 
       RowBox[{
       "u2x", ",", "u2y", ",", "u3x", ",", "u3y", ",", "u6x", ",", "u6y", ",",
         "u7x", ",", "u7y", ",", "u10x", ",", "u10y", ",", "u11x", ",", 
        "u11y", ",", "f1x", ",", "f1y", ",", "f4x", ",", "f4y", ",", "f5x", 
        ",", "f5y", ",", "f8x", ",", "f8y", ",", "f9x", ",", "f9y", ",", 
        "f12x", ",", "f12y"}], "}"}]}], "]"}]}], "\[IndentingNewLine]", 
   "VecIn", "\[IndentingNewLine]", 
   RowBox[{"Dimensions", "[", "VecIn", "]"}]}]}]], "Input", \
"PluginEmbeddedContent",
 CellGroupingRules->{GroupTogetherGrouping, 10000.}],

Cell[BoxData[
 RowBox[{"{", 
  RowBox[{"{", 
   RowBox[{
    RowBox[{"u2x", "\[Rule]", "1.2043269341790808`"}], ",", 
    RowBox[{"u2y", "\[Rule]", "0.4714106202926148`"}], ",", 
    RowBox[{"u3x", "\[Rule]", "2.178036339294035`"}], ",", 
    RowBox[{"u3y", "\[Rule]", "0.6178779072957218`"}], ",", 
    RowBox[{"u6x", "\[Rule]", "1.146469931345494`"}], ",", 
    RowBox[{"u6y", "\[Rule]", "0.07028440161321321`"}], ",", 
    RowBox[{"u7x", "\[Rule]", "2.4563650368616967`"}], ",", 
    RowBox[{"u7y", "\[Rule]", 
     RowBox[{"-", "0.010902849956788781`"}]}], ",", 
    RowBox[{"u10x", "\[Rule]", "1.6554477488184107`"}], ",", 
    RowBox[{"u10y", "\[Rule]", 
     RowBox[{"-", "0.5506141324093506`"}]}], ",", 
    RowBox[{"u11x", "\[Rule]", "2.4773356532845803`"}], ",", 
    RowBox[{"u11y", "\[Rule]", 
     RowBox[{"-", "0.6179363821314455`"}]}], ",", 
    RowBox[{"f1x", "\[Rule]", "8.972659526387224`*^6"}], ",", 
    RowBox[{"f1y", "\[Rule]", "2.3866271813889486`*^6"}], ",", 
    RowBox[{"f4x", "\[Rule]", 
     RowBox[{"-", "9.589476269822758`*^6"}]}], ",", 
    RowBox[{"f4y", "\[Rule]", "2.9746197370650377`*^6"}], ",", 
    RowBox[{"f5x", "\[Rule]", "1.8040934879913893`*^7"}], ",", 
    RowBox[{"f5y", "\[Rule]", "87594.92722199578`"}], ",", 
    RowBox[{"f8x", "\[Rule]", 
     RowBox[{"-", "1.6908921003903605`*^7"}]}], ",", 
    RowBox[{"f8y", "\[Rule]", 
     RowBox[{"-", "258216.63446468464`"}]}], ",", 
    RowBox[{"f9x", "\[Rule]", "8.14890351370029`*^6"}], ",", 
    RowBox[{"f9y", "\[Rule]", 
     RowBox[{"-", "2.5765516272988026`*^6"}]}], ",", 
    RowBox[{"f12x", "\[Rule]", 
     RowBox[{"-", "8.664100646275047`*^6"}]}], ",", 
    RowBox[{"f12y", "\[Rule]", 
     RowBox[{"-", "2.6140735839124992`*^6"}]}]}], "}"}], "}"}]], "Output", \
"PluginEmbeddedContent",
 CellGroupingRules->{GroupTogetherGrouping, 10000.}],

Cell[BoxData[
 RowBox[{"{", 
  RowBox[{
  "0", ",", "0", ",", "u2x", ",", "u2y", ",", "u3x", ",", "u3y", ",", "4.`", 
   ",", "0", ",", "0", ",", "0", ",", "u6x", ",", "u6y", ",", "u7x", ",", 
   "u7y", ",", "4.`", ",", "0", ",", "0", ",", "0", ",", "u10x", ",", "u10y", 
   ",", "u11x", ",", "u11y", ",", "4.`", ",", "0", ",", "f1x", ",", "f1y", 
   ",", "0", ",", "0", ",", "0", ",", "0", ",", "f4x", ",", "f4y", ",", "f5x",
    ",", "f5y", ",", "0", ",", "0", ",", "0", ",", "0", ",", "f8x", ",", 
   "f8y", ",", "f9x", ",", "f9y", ",", "0", ",", "0", ",", "0", ",", "0", ",",
    "f12x", ",", "f12y"}], "}"}]], "Output", "PluginEmbeddedContent",
 CellGroupingRules->{GroupTogetherGrouping, 10000.}],

Cell[BoxData[
 RowBox[{"{", "48", "}"}]], "Output", "PluginEmbeddedContent",
 CellGroupingRules->{GroupTogetherGrouping, 10000.}],

Cell[BoxData[{
 RowBox[{"Module", "[", 
  RowBox[{
   RowBox[{"{", 
    RowBox[{"j", ",", "i"}], "}"}], ",", 
   RowBox[{"(", "\[IndentingNewLine]", 
    RowBox[{
     RowBox[{"i", "=", "1"}], ";", "\[IndentingNewLine]", 
     RowBox[{"For", "[", 
      RowBox[{
       RowBox[{"j", "=", "1"}], ",", 
       RowBox[{"j", "<", 
        RowBox[{"(", " ", 
         RowBox[{"2", " ", "NNodos"}], ")"}]}], ",", 
       RowBox[{"j", "=", 
        RowBox[{"j", "+", "2"}]}], ",", 
       RowBox[{"(", "\[IndentingNewLine]", 
        RowBox[{
         RowBox[{
          RowBox[{"DeltaU", "[", 
           RowBox[{"[", 
            RowBox[{"i", ",", "1"}], "]"}], "]"}], "=", 
          RowBox[{"VecIn", "[", 
           RowBox[{"[", "j", "]"}], "]"}]}], ";", "\[IndentingNewLine]", 
         RowBox[{"i", "++"}], ";"}], "\[IndentingNewLine]", ")"}]}], "]"}], 
     ";", "\[IndentingNewLine]", 
     RowBox[{"i", "=", "1"}], ";", "\[IndentingNewLine]", 
     RowBox[{"For", "[", 
      RowBox[{
       RowBox[{"j", "=", "2"}], ",", 
       RowBox[{"j", "<", 
        RowBox[{"(", 
         RowBox[{
          RowBox[{"2", "  ", "NNodos"}], "+", "1"}], ")"}]}], ",", 
       RowBox[{"j", "=", 
        RowBox[{"j", "+", "2"}]}], ",", 
       RowBox[{"(", "\[IndentingNewLine]", 
        RowBox[{
         RowBox[{
          RowBox[{"DeltaU", "[", 
           RowBox[{"[", 
            RowBox[{"i", ",", "2"}], "]"}], "]"}], "=", 
          RowBox[{"VecIn", "[", 
           RowBox[{"[", "j", "]"}], "]"}]}], ";", "\[IndentingNewLine]", 
         RowBox[{"i", "++"}], ";"}], "\[IndentingNewLine]", ")"}]}], "]"}], 
     ";"}], "\[IndentingNewLine]", ")"}]}], 
  "]"}], "\[IndentingNewLine]", "DeltaU", "\[IndentingNewLine]", 
 RowBox[{"VecIn", "[", 
  RowBox[{"[", "3", "]"}], "]"}]}], "Input", "PluginEmbeddedContent",
 CellGroupingRules->{GroupTogetherGrouping, 10000.}],

Cell[BoxData[
 RowBox[{"{", 
  RowBox[{
   RowBox[{"{", 
    RowBox[{"0", ",", "0"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"u2x", ",", "u2y"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"u3x", ",", "u3y"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"4.`", ",", "0"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"0", ",", "0"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"u6x", ",", "u6y"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"u7x", ",", "u7y"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"4.`", ",", "0"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"0", ",", "0"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"u10x", ",", "u10y"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"u11x", ",", "u11y"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"4.`", ",", "0"}], "}"}]}], "}"}]], "Output", \
"PluginEmbeddedContent",
 CellGroupingRules->{GroupTogetherGrouping, 10000.}],

Cell[BoxData["u2x"], "Output", "PluginEmbeddedContent",
 CellGroupingRules->{GroupTogetherGrouping, 10000.}],

Cell[BoxData[{
 RowBox[{
  RowBox[{"\[Del]", "U"}], "=", 
  RowBox[{"DeltaU", "/.", 
   RowBox[{"First", "[", "Solucion", "]"}]}]}], "\[IndentingNewLine]", 
 RowBox[{"deformada", "=", 
  RowBox[{"p0", "+", 
   RowBox[{"\[Del]", "U"}]}]}], "\n", "p0", "\[IndentingNewLine]", 
 RowBox[{"MallaQ42", "[", 
  RowBox[{"4", ",", "3", ",", "p0"}], "]"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{"MallaDeform", "[", 
   RowBox[{"4", ",", "3", ",", "p0", ",", "deformada"}], "]"}], 
  "\[IndentingNewLine]", "\[IndentingNewLine]"}], "\n"}], "Input", \
"PluginEmbeddedContent",
 CellGroupingRules->{GroupTogetherGrouping, 10000.}],

Cell[BoxData[
 RowBox[{"{", 
  RowBox[{
   RowBox[{"{", 
    RowBox[{"0", ",", "0"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"1.2043269341790808`", ",", "0.4714106202926148`"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"2.178036339294035`", ",", "0.6178779072957218`"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"4.`", ",", "0"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"0", ",", "0"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"1.146469931345494`", ",", "0.07028440161321321`"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"2.4563650368616967`", ",", 
     RowBox[{"-", "0.010902849956788781`"}]}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"4.`", ",", "0"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"0", ",", "0"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"1.6554477488184107`", ",", 
     RowBox[{"-", "0.5506141324093506`"}]}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"2.4773356532845803`", ",", 
     RowBox[{"-", "0.6179363821314455`"}]}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"4.`", ",", "0"}], "}"}]}], "}"}]], "Output", \
"PluginEmbeddedContent",
 CellGroupingRules->{GroupTogetherGrouping, 10000.}],

Cell[BoxData[
 RowBox[{"{", 
  RowBox[{
   RowBox[{"{", 
    RowBox[{"0.`", ",", "0.`"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"4.204326934179081`", ",", "1.4714106202926147`"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"7.178036339294035`", ",", "1.1178779072957217`"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"14.`", ",", "0.`"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"0.`", ",", "5.`"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"4.146469931345494`", ",", "4.070284401613213`"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"8.456365036861698`", ",", "4.489097150043211`"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"14.`", ",", "5.`"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"0.`", ",", "9.`"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"5.6554477488184105`", ",", "7.9493858675906495`"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"8.47733565328458`", ",", "8.382063617868555`"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"14.`", ",", "9.`"}], "}"}]}], "}"}]], "Output", \
"PluginEmbeddedContent",
 CellGroupingRules->{GroupTogetherGrouping, 10000.}],

Cell[BoxData[
 RowBox[{"{", 
  RowBox[{
   RowBox[{"{", 
    RowBox[{"0.`", ",", "0.`"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"3.`", ",", "1.`"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"5", ",", "0.5`"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"10.`", ",", "0.`"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"0.`", ",", "5.`"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"3.`", ",", "4.`"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"6", ",", "4.5`"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"10.`", ",", "5.`"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"0.`", ",", "9.`"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"4.`", ",", "8.5`"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"6.`", ",", "9.`"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"10.`", ",", "9.`"}], "}"}]}], "}"}]], "Output", \
"PluginEmbeddedContent",
 CellGroupingRules->{GroupTogetherGrouping, 10000.}],

Cell[BoxData[
 GraphicsBox[
  GraphicsComplexBox[{{0., 0.}, {3., 1.}, {5, 0.5}, {10., 0.}, {0., 5.}, {3., 
    4.}, {6, 4.5}, {10., 5.}, {0., 9.}, {4., 8.5}, {6., 9.}, {10., 9.}}, {
    LineBox[{{1, 2, 6, 5, 1}, {2, 3, 7, 6, 2}, {3, 4, 8, 7, 3}, {5, 6, 10, 9, 
     5}, {6, 7, 11, 10, 6}, {7, 8, 12, 11, 7}}], 
    {RGBColor[0, 0, 1], PointSize[Large], 
     PointBox[{{0., 0.}, {3., 1.}, {5, 0.5}, {10., 0.}, {0., 5.}, {3., 4.}, {
       6, 4.5}, {10., 5.}, {0., 9.}, {4., 8.5}, {6., 9.}, {10., 
       9.}}]}}]]], "Output", "PluginEmbeddedContent",
 CellGroupingRules->{GroupTogetherGrouping, 10000.}],

Cell[BoxData[
 GraphicsBox[{
   GraphicsComplexBox[{{0., 0.}, {3., 1.}, {5, 0.5}, {10., 0.}, {0., 5.}, {3.,
      4.}, {6, 4.5}, {10., 5.}, {0., 9.}, {4., 8.5}, {6., 9.}, {10., 9.}}, {
     LineBox[{{1, 2, 6, 5, 1}, {2, 3, 7, 6, 2}, {3, 4, 8, 7, 3}, {5, 6, 10, 9,
       5}, {6, 7, 11, 10, 6}, {7, 8, 12, 11, 7}}], 
     {RGBColor[0, 0, 1], PointSize[Large], 
      PointBox[{{0., 0.}, {3., 1.}, {5, 0.5}, {10., 0.}, {0., 5.}, {3., 4.}, {
        6, 4.5}, {10., 5.}, {0., 9.}, {4., 8.5}, {6., 9.}, {10., 9.}}]}}], 
   GraphicsComplexBox[{{0., 0.}, {4.204326934179081, 1.4714106202926147`}, {
    7.178036339294035, 1.1178779072957217`}, {14., 0.}, {0., 5.}, {
    4.146469931345494, 4.070284401613213}, {8.456365036861698, 
    4.489097150043211}, {14., 5.}, {0., 9.}, {5.6554477488184105`, 
    7.9493858675906495`}, {8.47733565328458, 8.382063617868555}, {14., 
    9.}}, {LineBox[{{1, 2, 6, 5, 1}, {2, 3, 7, 6, 2}, {3, 4, 8, 7, 3}, {5, 6, 
      10, 9, 5}, {6, 7, 11, 10, 6}, {7, 8, 12, 11, 7}}], 
     {RGBColor[1, 0, 0], PointSize[Large], 
      PointBox[{{0., 0.}, {4.204326934179081, 1.4714106202926147`}, {
       7.178036339294035, 1.1178779072957217`}, {14., 0.}, {0., 5.}, {
       4.146469931345494, 4.070284401613213}, {8.456365036861698, 
       4.489097150043211}, {14., 5.}, {0., 9.}, {5.6554477488184105`, 
       7.9493858675906495`}, {8.47733565328458, 8.382063617868555}, {14., 
       9.}}]}}]}]], "Output", "PluginEmbeddedContent",
 CellGroupingRules->{GroupTogetherGrouping, 10000.}],

Cell["\<\
Si ahora aplicamos otro desplazamiento a la configuraci\[OAcute]n deformada:\
\>", "Text", "PluginEmbeddedContent",
 CellGroupingRules->{GroupTogetherGrouping, 10000.}],

Cell[BoxData[{
 RowBox[{
  RowBox[{"pref", "=", "p0"}], ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{"p0", "=", "deformada"}], ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{"GlobalMK", "[", "p0", "]"}], ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{"GlobalKAmp", "=", 
   RowBox[{"Array", "[", 
    RowBox[{
     RowBox[{"0", "&"}], ",", 
     RowBox[{"{", 
      RowBox[{
       RowBox[{"2", " ", "NNodos"}], ",", 
       RowBox[{"4", " ", "NNodos"}]}], "}"}]}], "]"}]}], 
  ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{
   RowBox[{"GlobalKAmp", "[", 
    RowBox[{"[", 
     RowBox[{"{", 
      RowBox[{"1", ",", "12"}], "}"}], "]"}], "]"}], "//", "MatrixForm"}], 
  ";"}], "\[IndentingNewLine]", 
 RowBox[{"For", "[", 
  RowBox[{
   RowBox[{"i", "=", "1"}], ",", 
   RowBox[{"i", "<", 
    RowBox[{"(", 
     RowBox[{
      RowBox[{"2", " ", "NNodos"}], "+", "1"}], ")"}]}], ",", 
   RowBox[{"i", "++"}], ",", 
   RowBox[{"(", "\[IndentingNewLine]", 
    RowBox[{"For", "[", 
     RowBox[{
      RowBox[{"j", "=", "1"}], ",", 
      RowBox[{"j", "<", 
       RowBox[{"(", 
        RowBox[{
         RowBox[{"2", " ", "NNodos"}], "+", "1"}], ")"}]}], ",", 
      RowBox[{"j", "++"}], ",", 
      RowBox[{"(", "\[IndentingNewLine]", 
       RowBox[{
        RowBox[{"GlobalKAmp", "[", 
         RowBox[{"[", 
          RowBox[{"i", ",", "j"}], "]"}], "]"}], "=", 
        RowBox[{"GlobalK", "[", 
         RowBox[{"[", 
          RowBox[{"i", ",", "j"}], "]"}], "]"}]}], "\[IndentingNewLine]", 
       ")"}]}], "]"}], "\[IndentingNewLine]", ")"}]}], 
  "]"}], "\[IndentingNewLine]", 
 RowBox[{"For", "[", 
  RowBox[{
   RowBox[{"i", "=", "1"}], ",", 
   RowBox[{"i", "<", 
    RowBox[{"(", 
     RowBox[{
      RowBox[{"2", " ", "NNodos"}], "+", "1"}], ")"}]}], ",", 
   RowBox[{"i", "++"}], ",", 
   RowBox[{"(", "\[IndentingNewLine]", 
    RowBox[{
     RowBox[{
      RowBox[{"GlobalKAmp", "[", 
       RowBox[{"[", 
        RowBox[{"i", ",", 
         RowBox[{"i", "+", 
          RowBox[{"2", " ", "NNodos"}]}]}], "]"}], "]"}], "=", "1"}], ";"}], 
    "\[IndentingNewLine]", ")"}]}], "]"}], "\n", 
 RowBox[{
  RowBox[{"GlobalKAmp", "//", "MatrixForm"}], ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{"Dimensions", "[", "GlobalKAmp", "]"}], 
  ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{"VecIn", "=", 
   RowBox[{"Array", "[", 
    RowBox[{
     RowBox[{"0", "&"}], ",", 
     RowBox[{"{", 
      RowBox[{"4", " ", "NNodos"}], "}"}]}], "]"}]}], 
  ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{"DeltaU", "=", 
   RowBox[{"Array", "[", 
    RowBox[{
     RowBox[{"0", "&"}], ",", 
     RowBox[{"{", " ", 
      RowBox[{"NNodos", ",", "2"}], "}"}]}], "]"}]}], 
  "\[IndentingNewLine]"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{
   RowBox[{"VecIn", "[", 
    RowBox[{"[", "7", "]"}], "]"}], "=", "4.0"}], ";", 
  RowBox[{
   RowBox[{"VecIn", "[", 
    RowBox[{"[", "15", "]"}], "]"}], "=", "4.0"}], ";", 
  RowBox[{
   RowBox[{"VecIn", "[", 
    RowBox[{"[", "23", "]"}], "]"}], "=", "4.0"}], 
  ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{
   RowBox[{
    RowBox[{"VecIn", "[", 
     RowBox[{"[", "3", "]"}], "]"}], "=", "u2x"}], ";", 
   RowBox[{
    RowBox[{"VecIn", "[", 
     RowBox[{"[", "4", "]"}], "]"}], "=", "u2y"}], ";", 
   RowBox[{
    RowBox[{"VecIn", "[", 
     RowBox[{"[", "5", "]"}], "]"}], "=", "u3x"}], ";", 
   RowBox[{
    RowBox[{"VecIn", "[", 
     RowBox[{"[", "6", "]"}], "]"}], "=", "u3y"}], ";", 
   RowBox[{
    RowBox[{"VecIn", "[", 
     RowBox[{"[", "11", "]"}], "]"}], "=", "u6x"}], ";", 
   RowBox[{
    RowBox[{"VecIn", "[", 
     RowBox[{"[", "12", "]"}], "]"}], "=", "u6y"}], ";", 
   RowBox[{
    RowBox[{"VecIn", "[", 
     RowBox[{"[", "13", "]"}], "]"}], "=", "u7x"}], ";", 
   RowBox[{
    RowBox[{"VecIn", "[", 
     RowBox[{"[", "14", "]"}], "]"}], "=", "u7y"}], ";", 
   RowBox[{
    RowBox[{"VecIn", "[", 
     RowBox[{"[", "19", "]"}], "]"}], "=", "u10x"}], ";", 
   RowBox[{
    RowBox[{"VecIn", "[", 
     RowBox[{"[", "20", "]"}], "]"}], "=", "u10y"}], ";", 
   RowBox[{
    RowBox[{"VecIn", "[", 
     RowBox[{"[", "21", "]"}], "]"}], "=", "u11x"}], ";", 
   RowBox[{
    RowBox[{"VecIn", "[", 
     RowBox[{"[", "22", "]"}], "]"}], "=", "u11y"}], ";"}], 
  "\[IndentingNewLine]"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{
   RowBox[{
    RowBox[{"VecIn", "[", 
     RowBox[{"[", "25", "]"}], "]"}], "=", "f1x"}], ";", 
   RowBox[{
    RowBox[{"VecIn", "[", 
     RowBox[{"[", "26", "]"}], "]"}], "=", "f1y"}], ";", 
   RowBox[{
    RowBox[{"VecIn", "[", 
     RowBox[{"[", "31", "]"}], "]"}], "=", "f4x"}], ";", 
   RowBox[{
    RowBox[{"VecIn", "[", 
     RowBox[{"[", "32", "]"}], "]"}], "=", "f4y"}], ";", 
   RowBox[{
    RowBox[{"VecIn", "[", 
     RowBox[{"[", "33", "]"}], "]"}], "=", "f5x"}], ";", 
   RowBox[{
    RowBox[{"VecIn", "[", 
     RowBox[{"[", "34", "]"}], "]"}], "=", "f5y"}], ";", 
   RowBox[{
    RowBox[{"VecIn", "[", 
     RowBox[{"[", "39", "]"}], "]"}], "=", "f8x"}], ";", 
   RowBox[{
    RowBox[{"VecIn", "[", 
     RowBox[{"[", "40", "]"}], "]"}], "=", "f8y"}], ";", 
   RowBox[{
    RowBox[{"VecIn", "[", 
     RowBox[{"[", "41", "]"}], "]"}], "=", "f9x"}], ";", 
   RowBox[{
    RowBox[{"VecIn", "[", 
     RowBox[{"[", "42", "]"}], "]"}], "=", "f9y"}], ";", 
   RowBox[{
    RowBox[{"VecIn", "[", 
     RowBox[{"[", "47", "]"}], "]"}], "=", "f12x"}], ";", 
   RowBox[{
    RowBox[{"VecIn", "[", 
     RowBox[{"[", "48", "]"}], "]"}], "=", "f12y"}], ";"}], 
  "\[IndentingNewLine]"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{"Solucion", "=", 
   RowBox[{"GlobalKAmp", ".", "VecIn"}]}], ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{"Dimensions", "[", "Solucion", "]"}], ";"}], "\[IndentingNewLine]", 

 RowBox[{
  RowBox[{"Solucion", "=", 
   RowBox[{"NSolve", "[", 
    RowBox[{
     RowBox[{"{", 
      RowBox[{
       RowBox[{
        RowBox[{"Solucion", "[", 
         RowBox[{"[", "2", "]"}], "]"}], "\[Equal]", "0"}], " ", ",", " ", 
       RowBox[{
        RowBox[{"Solucion", "[", 
         RowBox[{"[", "1", "]"}], "]"}], "\[Equal]", "0"}], " ", ",", 
       RowBox[{
        RowBox[{"Solucion", "[", 
         RowBox[{"[", "3", "]"}], "]"}], "\[Equal]", "0"}], ",", 
       RowBox[{
        RowBox[{"Solucion", "[", 
         RowBox[{"[", "4", "]"}], "]"}], "\[Equal]", "0"}], ",", 
       RowBox[{
        RowBox[{"Solucion", "[", 
         RowBox[{"[", "5", "]"}], "]"}], "\[Equal]", "0"}], 
       "\[IndentingNewLine]", ",", 
       RowBox[{
        RowBox[{"Solucion", "[", 
         RowBox[{"[", "6", "]"}], "]"}], "\[Equal]", "0"}], ",", 
       RowBox[{
        RowBox[{"Solucion", "[", 
         RowBox[{"[", "7", "]"}], "]"}], "\[Equal]", "0"}], " ", ",", 
       RowBox[{
        RowBox[{"Solucion", "[", 
         RowBox[{"[", "9", "]"}], "]"}], "\[Equal]", "0"}], " ", ",", 
       RowBox[{
        RowBox[{"Solucion", "[", 
         RowBox[{"[", "8", "]"}], "]"}], "\[Equal]", "0"}], " ", ",", 
       RowBox[{
        RowBox[{"Solucion", "[", 
         RowBox[{"[", "10", "]"}], "]"}], "\[Equal]", "0"}], 
       "\[IndentingNewLine]", ",", 
       RowBox[{
        RowBox[{"Solucion", "[", 
         RowBox[{"[", "11", "]"}], "]"}], "\[Equal]", "0"}], ",", 
       RowBox[{
        RowBox[{"Solucion", "[", 
         RowBox[{"[", "12", "]"}], "]"}], "\[Equal]", "0"}], ",", 
       RowBox[{
        RowBox[{"Solucion", "[", 
         RowBox[{"[", "13", "]"}], "]"}], "\[Equal]", "0"}], ",", 
       RowBox[{
        RowBox[{"Solucion", "[", 
         RowBox[{"[", "14", "]"}], "]"}], "\[Equal]", "0"}], " ", ",", 
       RowBox[{
        RowBox[{"Solucion", "[", 
         RowBox[{"[", "15", "]"}], "]"}], "\[Equal]", "0"}], " ", ",", 
       RowBox[{
        RowBox[{"Solucion", "[", 
         RowBox[{"[", "16", "]"}], "]"}], "\[Equal]", "0"}], "  ", ",", 
       RowBox[{
        RowBox[{"Solucion", "[", 
         RowBox[{"[", "17", "]"}], "]"}], "\[Equal]", "0"}], " ", ",", 
       RowBox[{
        RowBox[{"Solucion", "[", 
         RowBox[{"[", "18", "]"}], "]"}], "\[Equal]", "0"}], " ", ",", 
       RowBox[{
        RowBox[{"Solucion", "[", 
         RowBox[{"[", "19", "]"}], "]"}], "\[Equal]", "0"}], ",", 
       RowBox[{
        RowBox[{"Solucion", "[", 
         RowBox[{"[", "20", "]"}], "]"}], "\[Equal]", "0"}], ",", 
       RowBox[{
        RowBox[{"Solucion", "[", 
         RowBox[{"[", "21", "]"}], "]"}], "\[Equal]", "0"}], ",", 
       RowBox[{
        RowBox[{"Solucion", "[", 
         RowBox[{"[", "22", "]"}], "]"}], "\[Equal]", "0"}], " ", ",", 
       RowBox[{
        RowBox[{"Solucion", "[", 
         RowBox[{"[", "23", "]"}], "]"}], "\[Equal]", "0"}], " ", ",", 
       RowBox[{
        RowBox[{"Solucion", "[", 
         RowBox[{"[", "24", "]"}], "]"}], "\[Equal]", "0"}]}], " ", " ", " ", 
      " ", "}"}], ",", 
     RowBox[{"{", 
      RowBox[{
      "u2x", ",", "u2y", ",", "u3x", ",", "u3y", ",", "u6x", ",", "u6y", ",", 
       "u7x", ",", "u7y", ",", "u10x", ",", "u10y", ",", "u11x", ",", "u11y", 
       ",", "f1x", ",", "f1y", ",", "f4x", ",", "f4y", ",", "f5x", ",", "f5y",
        ",", "f8x", ",", "f8y", ",", "f9x", ",", "f9y", ",", "f12x", ",", 
       "f12y"}], "}"}]}], "]"}]}], 
  "\[IndentingNewLine]"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{"Module", "[", 
   RowBox[{
    RowBox[{"{", 
     RowBox[{"j", ",", "i"}], "}"}], ",", 
    RowBox[{"(", "\[IndentingNewLine]", 
     RowBox[{
      RowBox[{"i", "=", "1"}], ";", "\[IndentingNewLine]", 
      RowBox[{"For", "[", 
       RowBox[{
        RowBox[{"j", "=", "1"}], ",", 
        RowBox[{"j", "<", 
         RowBox[{"(", " ", 
          RowBox[{"2", " ", "NNodos"}], ")"}]}], ",", 
        RowBox[{"j", "=", 
         RowBox[{"j", "+", "2"}]}], ",", 
        RowBox[{"(", "\[IndentingNewLine]", 
         RowBox[{
          RowBox[{
           RowBox[{"DeltaU", "[", 
            RowBox[{"[", 
             RowBox[{"i", ",", "1"}], "]"}], "]"}], "=", 
           RowBox[{"VecIn", "[", 
            RowBox[{"[", "j", "]"}], "]"}]}], ";", "\[IndentingNewLine]", 
          RowBox[{"i", "++"}], ";"}], "\[IndentingNewLine]", ")"}]}], "]"}], 
      ";", "\[IndentingNewLine]", 
      RowBox[{"i", "=", "1"}], ";", "\[IndentingNewLine]", 
      RowBox[{"For", "[", 
       RowBox[{
        RowBox[{"j", "=", "2"}], ",", 
        RowBox[{"j", "<", 
         RowBox[{"(", 
          RowBox[{
           RowBox[{"2", "  ", "NNodos"}], "+", "1"}], ")"}]}], ",", 
        RowBox[{"j", "=", 
         RowBox[{"j", "+", "2"}]}], ",", 
        RowBox[{"(", "\[IndentingNewLine]", 
         RowBox[{
          RowBox[{
           RowBox[{"DeltaU", "[", 
            RowBox[{"[", 
             RowBox[{"i", ",", "2"}], "]"}], "]"}], "=", 
           RowBox[{"VecIn", "[", 
            RowBox[{"[", "j", "]"}], "]"}]}], ";", "\[IndentingNewLine]", 
          RowBox[{"i", "++"}], ";"}], "\[IndentingNewLine]", ")"}]}], "]"}], 
      ";"}], "\[IndentingNewLine]", ")"}]}], "]"}], 
  "\[IndentingNewLine]"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{"\[Del]", "U"}], "=", 
  RowBox[{"DeltaU", "/.", 
   RowBox[{"First", "[", "Solucion", "]"}]}]}], "\[IndentingNewLine]", 
 RowBox[{"deformada", "=", 
  RowBox[{"p0", "+", 
   RowBox[{"\[Del]", "U"}]}]}], "\n", "p0", "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{"MallaDeform", "[", 
   RowBox[{"4", ",", "3", ",", "pref", ",", "deformada"}], "]"}], 
  "\[IndentingNewLine]", "\[IndentingNewLine]", 
  "\[IndentingNewLine]"}], "\[IndentingNewLine]"}], "Input", "PluginEmbeddedContent",
 CellGroupingRules->{GroupTogetherGrouping, 10000.}],

Cell[BoxData[
 RowBox[{"{", 
  RowBox[{
   RowBox[{"{", 
    RowBox[{"0", ",", "0"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"0", ",", "0"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"0", ",", "0"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"0", ",", "0"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"0", ",", "0"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"0", ",", "0"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"0", ",", "0"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"0", ",", "0"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"0", ",", "0"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"0", ",", "0"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"0", ",", "0"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"0", ",", "0"}], "}"}]}], "}"}]], "Output", \
"PluginEmbeddedContent",
 CellGroupingRules->{GroupTogetherGrouping, 10000.}],

Cell[BoxData[
 RowBox[{"{", 
  RowBox[{"{", 
   RowBox[{
    RowBox[{"u2x", "\[Rule]", "1.1036221401321642`"}], ",", 
    RowBox[{"u2y", "\[Rule]", "0.32609285225053847`"}], ",", 
    RowBox[{"u3x", "\[Rule]", "2.2593377512733888`"}], ",", 
    RowBox[{"u3y", "\[Rule]", "0.4086979367749504`"}], ",", 
    RowBox[{"u6x", "\[Rule]", "1.1469830201704803`"}], ",", 
    RowBox[{"u6y", "\[Rule]", "0.04459246405661588`"}], ",", 
    RowBox[{"u7x", "\[Rule]", "2.492080946148787`"}], ",", 
    RowBox[{"u7y", "\[Rule]", 
     RowBox[{"-", "0.01782971637822274`"}]}], ",", 
    RowBox[{"u10x", "\[Rule]", "1.5984309587632004`"}], ",", 
    RowBox[{"u10y", "\[Rule]", 
     RowBox[{"-", "0.4049464559199114`"}]}], ",", 
    RowBox[{"u11x", "\[Rule]", "2.5484131340954685`"}], ",", 
    RowBox[{"u11y", "\[Rule]", 
     RowBox[{"-", "0.4291456512068984`"}]}], ",", 
    RowBox[{"f1x", "\[Rule]", "5.442330146736493`*^6"}], ",", 
    RowBox[{"f1y", "\[Rule]", "1.9291500280603964`*^6"}], ",", 
    RowBox[{"f4x", "\[Rule]", 
     RowBox[{"-", "5.886237893686771`*^6"}]}], ",", 
    RowBox[{"f4y", "\[Rule]", "2.4080538323283526`*^6"}], ",", 
    RowBox[{"f5x", "\[Rule]", "1.250057313119493`*^7"}], ",", 
    RowBox[{"f5y", "\[Rule]", "132117.65439890785`"}], ",", 
    RowBox[{"f8x", "\[Rule]", 
     RowBox[{"-", "1.1723192298131915`*^7"}]}], ",", 
    RowBox[{"f8y", "\[Rule]", 
     RowBox[{"-", "263497.2991260508`"}]}], ",", 
    RowBox[{"f9x", "\[Rule]", "4.849235559713391`*^6"}], ",", 
    RowBox[{"f9y", "\[Rule]", 
     RowBox[{"-", "2.124528138909327`*^6"}]}], ",", 
    RowBox[{"f12x", "\[Rule]", 
     RowBox[{"-", "5.182708645826147`*^6"}]}], ",", 
    RowBox[{"f12y", "\[Rule]", 
     RowBox[{"-", "2.0812960767522727`*^6"}]}]}], "}"}], "}"}]], "Output", \
"PluginEmbeddedContent",
 CellGroupingRules->{GroupTogetherGrouping, 10000.}],

Cell[BoxData[
 RowBox[{"{", 
  RowBox[{
   RowBox[{"{", 
    RowBox[{"0", ",", "0"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"1.1036221401321642`", ",", "0.32609285225053847`"}], "}"}], ",", 
   
   RowBox[{"{", 
    RowBox[{"2.2593377512733888`", ",", "0.4086979367749504`"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"4.`", ",", "0"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"0", ",", "0"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"1.1469830201704803`", ",", "0.04459246405661588`"}], "}"}], ",", 
   
   RowBox[{"{", 
    RowBox[{"2.492080946148787`", ",", 
     RowBox[{"-", "0.01782971637822274`"}]}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"4.`", ",", "0"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"0", ",", "0"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"1.5984309587632004`", ",", 
     RowBox[{"-", "0.4049464559199114`"}]}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"2.5484131340954685`", ",", 
     RowBox[{"-", "0.4291456512068984`"}]}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"4.`", ",", "0"}], "}"}]}], "}"}]], "Output", \
"PluginEmbeddedContent",
 CellGroupingRules->{GroupTogetherGrouping, 10000.}],

Cell[BoxData[
 RowBox[{"{", 
  RowBox[{
   RowBox[{"{", 
    RowBox[{"0.`", ",", "0.`"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"5.307949074311245`", ",", "1.7975034725431531`"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"9.437374090567424`", ",", "1.526575844070672`"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"18.`", ",", "0.`"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"0.`", ",", "5.`"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"5.293452951515974`", ",", "4.114876865669829`"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"10.948445983010485`", ",", "4.4712674336649885`"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"18.`", ",", "5.`"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"0.`", ",", "9.`"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"7.253878707581611`", ",", "7.544439411670738`"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"11.02574878738005`", ",", "7.952917966661657`"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"18.`", ",", "9.`"}], "}"}]}], "}"}]], "Output", \
"PluginEmbeddedContent",
 CellGroupingRules->{GroupTogetherGrouping, 10000.}],

Cell[BoxData[
 RowBox[{"{", 
  RowBox[{
   RowBox[{"{", 
    RowBox[{"0.`", ",", "0.`"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"4.204326934179081`", ",", "1.4714106202926147`"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"7.178036339294035`", ",", "1.1178779072957217`"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"14.`", ",", "0.`"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"0.`", ",", "5.`"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"4.146469931345494`", ",", "4.070284401613213`"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"8.456365036861698`", ",", "4.489097150043211`"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"14.`", ",", "5.`"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"0.`", ",", "9.`"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"5.6554477488184105`", ",", "7.9493858675906495`"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"8.47733565328458`", ",", "8.382063617868555`"}], "}"}], ",", 
   RowBox[{"{", 
    RowBox[{"14.`", ",", "9.`"}], "}"}]}], "}"}]], "Output", \
"PluginEmbeddedContent",
 CellGroupingRules->{GroupTogetherGrouping, 10000.}],

Cell[BoxData[
 GraphicsBox[{
   GraphicsComplexBox[{{0., 0.}, {3., 1.}, {5, 0.5}, {10., 0.}, {0., 5.}, {3.,
      4.}, {6, 4.5}, {10., 5.}, {0., 9.}, {4., 8.5}, {6., 9.}, {10., 9.}}, {
     LineBox[{{1, 2, 6, 5, 1}, {2, 3, 7, 6, 2}, {3, 4, 8, 7, 3}, {5, 6, 10, 9,
       5}, {6, 7, 11, 10, 6}, {7, 8, 12, 11, 7}}], 
     {RGBColor[0, 0, 1], PointSize[Large], 
      PointBox[{{0., 0.}, {3., 1.}, {5, 0.5}, {10., 0.}, {0., 5.}, {3., 4.}, {
        6, 4.5}, {10., 5.}, {0., 9.}, {4., 8.5}, {6., 9.}, {10., 9.}}]}}], 
   GraphicsComplexBox[{{0., 0.}, {5.307949074311245, 1.7975034725431531`}, {
    9.437374090567424, 1.526575844070672}, {18., 0.}, {0., 5.}, {
    5.293452951515974, 4.114876865669829}, {10.948445983010485`, 
    4.4712674336649885`}, {18., 5.}, {0., 9.}, {7.253878707581611, 
    7.544439411670738}, {11.02574878738005, 7.952917966661657}, {18., 
    9.}}, {LineBox[{{1, 2, 6, 5, 1}, {2, 3, 7, 6, 2}, {3, 4, 8, 7, 3}, {5, 6, 
      10, 9, 5}, {6, 7, 11, 10, 6}, {7, 8, 12, 11, 7}}], 
     {RGBColor[1, 0, 0], PointSize[Large], 
      PointBox[{{0., 0.}, {5.307949074311245, 1.7975034725431531`}, {
       9.437374090567424, 1.526575844070672}, {18., 0.}, {0., 5.}, {
       5.293452951515974, 4.114876865669829}, {10.948445983010485`, 
       4.4712674336649885`}, {18., 5.}, {0., 9.}, {7.253878707581611, 
       7.544439411670738}, {11.02574878738005, 7.952917966661657}, {18., 
       9.}}]}}]}]], "Output", "PluginEmbeddedContent",
 CellGroupingRules->{GroupTogetherGrouping, 10000.}]
}, Open  ]],

Cell[BoxData[{
 RowBox[{
  RowBox[{"p0", "=", 
   RowBox[{"{", 
    RowBox[{
     RowBox[{"{", 
      RowBox[{"0.", ",", "0."}], "}"}], ",", 
     RowBox[{"{", 
      RowBox[{"3.", ",", "1."}], "}"}], ",", 
     RowBox[{"{", 
      RowBox[{"10.", ",", "0."}], "}"}], ",", 
     RowBox[{"{", 
      RowBox[{"0.", ",", "5."}], "}"}], ",", 
     RowBox[{"{", 
      RowBox[{"3.", ",", "4."}], "}"}], ",", 
     RowBox[{"{", 
      RowBox[{"10.", ",", "5."}], "}"}]}], "}"}]}], 
  ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{"p1", "=", 
   RowBox[{"{", 
    RowBox[{
     RowBox[{"{", 
      RowBox[{"0.", ",", "0."}], "}"}], ",", 
     RowBox[{"{", 
      RowBox[{"3.", ",", "0."}], "}"}], ",", 
     RowBox[{"{", 
      RowBox[{"10.", ",", "0."}], "}"}], ",", "\[IndentingNewLine]", 
     RowBox[{"{", 
      RowBox[{"0.", ",", "5."}], "}"}], ",", 
     RowBox[{"{", 
      RowBox[{"3.", ",", "4."}], "}"}], ",", 
     RowBox[{"{", 
      RowBox[{"10.", ",", "5."}], "}"}], ",", "\[IndentingNewLine]", 
     RowBox[{"{", 
      RowBox[{"0.0", ",", "8.0"}], "}"}], ",", 
     RowBox[{"{", 
      RowBox[{"4.0", ",", "7.0"}], "}"}], ",", 
     RowBox[{"{", 
      RowBox[{"10.0", ",", "9.0"}], "}"}], ",", "\[IndentingNewLine]", 
     RowBox[{"{", 
      RowBox[{"0.0", ",", "12.0"}], "}"}], ",", 
     RowBox[{"{", 
      RowBox[{"5.0", ",", "12.0"}], "}"}], ",", 
     RowBox[{"{", 
      RowBox[{"10.0", ",", "12.0"}], "}"}]}], "}"}]}], 
  ";"}], "\[IndentingNewLine]", 
 RowBox[{"MallaQ42", "[", 
  RowBox[{"3", ",", "2", ",", "p0"}], "]"}], "\n"}], "Input", \
"PluginEmbeddedContent"]
}, Open  ]]
},
WindowSize->{2017, 9219},
Visible->True,
AuthoredSize->{2017, 9219},
ScrollingOptions->{"HorizontalScrollRange"->Fit,
"VerticalScrollRange"->Fit},
ShowCellBracket->False,
Deployed->True,
CellContext->Notebook,
TrackCellChangeTimes->False,
FrontEndVersion->"9.0 for Mac OS X x86 (32-bit, 64-bit Kernel) (November 20, \
2012)",
StyleDefinitions->"Default.nb"
]
(* End of Notebook Content *)

(* Internal cache information *)
(*CellTagsOutline
CellTagsIndex->{}
*)
(*CellTagsIndex
CellTagsIndex->{}
*)
(*NotebookFileOutline
Notebook[{
Cell[CellGroupData[{
Cell[1485, 35, 95, 2, 38, "Section"],
Cell[CellGroupData[{
Cell[1605, 41, 130, 3, 49, "Subsection"],
Cell[1738, 46, 76, 0, 16, "Text"],
Cell[1817, 48, 469, 12, 65, "Input"],
Cell[2289, 62, 64, 0, 16, "Text"],
Cell[CellGroupData[{
Cell[2378, 66, 510, 16, 13, "Input"],
Cell[2891, 84, 626, 16, 232, "Output"],
Cell[3520, 102, 680, 18, 227, "Output"]
}, Open  ]],
Cell[4215, 123, 112, 2, 17, "Text"],
Cell[4330, 127, 944, 33, 65, "Input"],
Cell[5277, 162, 75, 0, 16, "Text"],
Cell[CellGroupData[{
Cell[5377, 166, 658, 21, 31, "Input"],
Cell[6038, 189, 69444, 1135, 286, "Output"],
Cell[75485, 1326, 244516, 3929, 286, "Output"]
}, Open  ]],
Cell[320016, 5258, 182, 3, 17, "Text"],
Cell[CellGroupData[{
Cell[320223, 5265, 500, 14, 31, "Input"],
Cell[320726, 5281, 179, 4, 13, "Output"],
Cell[320908, 5287, 178, 4, 13, "Output"],
Cell[321089, 5293, 179, 4, 13, "Output"],
Cell[321271, 5299, 178, 4, 13, "Output"],
Cell[321452, 5305, 179, 4, 13, "Output"],
Cell[321634, 5311, 178, 4, 13, "Output"],
Cell[321815, 5317, 179, 4, 13, "Output"],
Cell[321997, 5323, 325, 8, 13, "Output"],
Cell[322325, 5333, 325, 8, 13, "Output"],
Cell[322653, 5343, 179, 4, 13, "Output"]
}, Open  ]],
Cell[322847, 5350, 107, 2, 16, "Text"],
Cell[CellGroupData[{
Cell[322979, 5356, 572, 18, 13, "Input"],
Cell[323554, 5376, 666, 19, 32, "Output"],
Cell[324223, 5397, 666, 19, 32, "Output"],
Cell[324892, 5418, 666, 19, 32, "Output"],
Cell[325561, 5439, 838, 23, 32, "Output"],
Cell[326402, 5464, 838, 23, 32, "Output"]
}, Open  ]],
Cell[327255, 5490, 88, 0, 16, "Text"],
Cell[CellGroupData[{
Cell[327368, 5494, 109, 2, 13, "Input"],
Cell[327480, 5498, 107, 2, 13, "Output"],
Cell[327590, 5502, 107, 2, 13, "Output"],
Cell[327700, 5506, 107, 2, 13, "Output"],
Cell[327810, 5510, 183, 4, 13, "Output"],
Cell[327996, 5516, 183, 4, 13, "Output"]
}, Open  ]],
Cell[328194, 5523, 100, 2, 16, "Text"],
Cell[CellGroupData[{
Cell[328319, 5529, 176, 5, 13, "Input"],
Cell[328498, 5536, 670, 19, 32, "Output"],
Cell[329171, 5557, 670, 19, 32, "Output"],
Cell[329844, 5578, 670, 19, 32, "Output"],
Cell[330517, 5599, 1197, 33, 54, "Output"],
Cell[331717, 5634, 1197, 33, 54, "Output"]
}, Open  ]],
Cell[332929, 5670, 91, 0, 16, "Text"],
Cell[CellGroupData[{
Cell[333045, 5674, 578, 17, 65, "Input"],
Cell[333626, 5693, 106, 2, 13, "Output"],
Cell[333735, 5697, 105, 2, 13, "Output"],
Cell[333843, 5701, 105, 2, 13, "Output"],
Cell[333951, 5705, 106, 2, 13, "Output"],
Cell[334060, 5709, 106, 2, 13, "Output"],
Cell[334169, 5713, 105, 2, 13, "Output"],
Cell[334277, 5717, 105, 2, 13, "Output"],
Cell[334385, 5721, 106, 2, 13, "Output"],
Cell[334494, 5725, 106, 2, 13, "Output"],
Cell[334603, 5729, 105, 2, 13, "Output"],
Cell[334711, 5733, 105, 2, 13, "Output"],
Cell[334819, 5737, 106, 2, 13, "Output"],
Cell[334928, 5741, 285, 7, 35, "Output"],
Cell[335216, 5750, 103, 2, 13, "Output"],
Cell[335322, 5754, 283, 8, 35, "Output"],
Cell[335608, 5764, 208, 5, 35, "Output"],
Cell[335819, 5771, 208, 5, 35, "Output"],
Cell[336030, 5778, 283, 8, 35, "Output"],
Cell[336316, 5788, 103, 2, 13, "Output"],
Cell[336422, 5792, 285, 7, 35, "Output"]
}, Open  ]],
Cell[336722, 5802, 148, 3, 16, "Text"],
Cell[CellGroupData[{
Cell[336895, 5809, 539, 18, 31, "Input"],
Cell[337437, 5829, 249, 8, 13, "Output"],
Cell[337689, 5839, 250, 8, 13, "Output"],
Cell[337942, 5849, 249, 8, 13, "Output"],
Cell[338194, 5859, 250, 8, 13, "Output"],
Cell[338447, 5869, 249, 8, 13, "Output"],
Cell[338699, 5879, 250, 8, 13, "Output"],
Cell[338952, 5889, 869, 29, 35, "Output"],
Cell[339824, 5920, 308, 9, 35, "Output"],
Cell[340135, 5931, 308, 9, 35, "Output"],
Cell[340446, 5942, 869, 29, 35, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[341352, 5976, 539, 18, 31, "Input"],
Cell[341894, 5996, 249, 8, 13, "Output"],
Cell[342146, 6006, 229, 7, 13, "Output"],
Cell[342378, 6015, 249, 8, 13, "Output"],
Cell[342630, 6025, 229, 7, 13, "Output"],
Cell[342862, 6034, 249, 8, 13, "Output"],
Cell[343114, 6044, 229, 7, 13, "Output"],
Cell[343346, 6053, 847, 28, 35, "Output"],
Cell[344196, 6083, 287, 8, 35, "Output"],
Cell[344486, 6093, 287, 8, 35, "Output"],
Cell[344776, 6103, 847, 28, 35, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[345660, 6136, 539, 18, 31, "Input"],
Cell[346202, 6156, 226, 6, 13, "Output"],
Cell[346431, 6164, 229, 7, 13, "Output"],
Cell[346663, 6173, 226, 6, 13, "Output"],
Cell[346892, 6181, 229, 7, 13, "Output"],
Cell[347124, 6190, 226, 6, 13, "Output"],
Cell[347353, 6198, 229, 7, 13, "Output"],
Cell[347585, 6207, 825, 27, 35, "Output"],
Cell[348413, 6236, 287, 8, 35, "Output"],
Cell[348703, 6246, 287, 8, 35, "Output"],
Cell[348993, 6256, 825, 27, 35, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[349855, 6288, 539, 18, 31, "Input"],
Cell[350397, 6308, 226, 6, 13, "Output"],
Cell[350626, 6316, 250, 8, 13, "Output"],
Cell[350879, 6326, 226, 6, 13, "Output"],
Cell[351108, 6334, 250, 8, 13, "Output"],
Cell[351361, 6344, 226, 6, 13, "Output"],
Cell[351590, 6352, 250, 8, 13, "Output"],
Cell[351843, 6362, 847, 28, 35, "Output"],
Cell[352693, 6392, 308, 9, 35, "Output"],
Cell[353004, 6403, 308, 9, 35, "Output"],
Cell[353315, 6414, 847, 28, 35, "Output"]
}, Open  ]],
Cell[354177, 6445, 71, 0, 16, "Text"],
Cell[CellGroupData[{
Cell[354273, 6449, 704, 18, 48, "Input"],
Cell[354980, 6469, 3271, 92, 48, "Output"],
Cell[358254, 6563, 3271, 92, 48, "Output"],
Cell[361528, 6657, 3271, 92, 48, "Output"],
Cell[364802, 6751, 9519, 272, 82, "Output"],
Cell[374324, 7025, 9519, 272, 82, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[383880, 7302, 824, 26, 48, "Input"],
Cell[384707, 7330, 803, 20, 54, "Output"],
Cell[385513, 7352, 108, 2, 17, "Output"],
Cell[385624, 7356, 803, 20, 54, "Output"],
Cell[386430, 7378, 108, 2, 17, "Output"],
Cell[386541, 7382, 803, 20, 54, "Output"],
Cell[387347, 7404, 108, 2, 17, "Output"],
Cell[387458, 7408, 803, 20, 54, "Output"],
Cell[388264, 7430, 108, 2, 17, "Output"],
Cell[388375, 7434, 108, 2, 17, "Output"],
Cell[388486, 7438, 803, 20, 54, "Output"]
}, Open  ]],
Cell[389304, 7461, 113, 2, 16, "Text"],
Cell[389420, 7465, 459, 14, 31, "Input"],
Cell[389882, 7481, 71, 0, 16, "Text"],
Cell[CellGroupData[{
Cell[389978, 7485, 1157, 34, 31, "Input"],
Cell[391138, 7521, 3002, 65, 150, "Output"],
Cell[394143, 7588, 3002, 65, 150, "Output"],
Cell[397148, 7655, 3002, 65, 150, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[400187, 7725, 1717, 54, 48, "Input"],
Cell[401907, 7781, 1835, 44, 150, "Output"],
Cell[403745, 7827, 1835, 44, 150, "Output"],
Cell[405583, 7873, 1835, 44, 150, "Output"],
Cell[407421, 7919, 1849, 44, 150, "Output"],
Cell[409273, 7965, 1848, 45, 150, "Output"]
}, Open  ]],
Cell[CellGroupData[{
Cell[411158, 8015, 1637, 52, 65, "Input"],
Cell[412798, 8069, 666, 16, 13, "Output"]
}, Open  ]]
}, Open  ]]
}, Open  ]],
Cell[CellGroupData[{
Cell[413525, 8092, 93, 2, 38, "Section"],
Cell[413621, 8096, 93, 0, 16, "Text"],
Cell[413717, 8098, 372, 10, 48, "Input"],
Cell[CellGroupData[{
Cell[414114, 8112, 186, 4, 13, "Input",
 CellGroupingRules->{GroupTogetherGrouping, 10000.}],
Cell[414303, 8118, 800, 13, 187, "Output",
 CellGroupingRules->{GroupTogetherGrouping, 10000.}],
Cell[415106, 8133, 157, 3, 16, "Text",
 CellGroupingRules->{GroupTogetherGrouping, 10000.}],
Cell[415266, 8138, 1391, 42, 82, "Input",
 CellGroupingRules->{GroupTogetherGrouping, 10000.}],
Cell[416660, 8182, 2546, 82, 235, "Input",
 CellGroupingRules->{GroupTogetherGrouping, 10000.}],
Cell[419209, 8266, 841, 28, 13, "Output",
 CellGroupingRules->{GroupTogetherGrouping, 10000.}],
Cell[420053, 8296, 188, 4, 16, "Text",
 CellGroupingRules->{GroupTogetherGrouping, 10000.}],
Cell[420244, 8302, 1563, 50, 65, "Input",
 CellGroupingRules->{GroupTogetherGrouping, 10000.}],
Cell[421810, 8354, 127, 1, 16, "Text",
 CellGroupingRules->{GroupTogetherGrouping, 10000.}],
Cell[421940, 8357, 1244, 39, 48, "Input",
 CellGroupingRules->{GroupTogetherGrouping, 10000.}],
Cell[423187, 8398, 4005, 99, 167, "Input",
 CellGroupingRules->{GroupTogetherGrouping, 10000.}],
Cell[427195, 8499, 1841, 38, 62, "Output",
 CellGroupingRules->{GroupTogetherGrouping, 10000.}],
Cell[429039, 8539, 702, 11, 31, "Output",
 CellGroupingRules->{GroupTogetherGrouping, 10000.}],
Cell[429744, 8552, 129, 2, 13, "Output",
 CellGroupingRules->{GroupTogetherGrouping, 10000.}],
Cell[429876, 8556, 1863, 49, 235, "Input",
 CellGroupingRules->{GroupTogetherGrouping, 10000.}],
Cell[431742, 8607, 875, 28, 13, "Output",
 CellGroupingRules->{GroupTogetherGrouping, 10000.}],
Cell[432620, 8637, 108, 1, 13, "Output",
 CellGroupingRules->{GroupTogetherGrouping, 10000.}],
Cell[432731, 8640, 623, 15, 133, "Input",
 CellGroupingRules->{GroupTogetherGrouping, 10000.}],
Cell[433357, 8657, 1127, 31, 31, "Output",
 CellGroupingRules->{GroupTogetherGrouping, 10000.}],
Cell[434487, 8690, 1075, 28, 31, "Output",
 CellGroupingRules->{GroupTogetherGrouping, 10000.}],
Cell[435565, 8720, 891, 28, 13, "Output",
 CellGroupingRules->{GroupTogetherGrouping, 10000.}],
Cell[436459, 8750, 603, 10, 325, "Output",
 CellGroupingRules->{GroupTogetherGrouping, 10000.}],
Cell[437065, 8762, 1506, 23, 236, "Output",
 CellGroupingRules->{GroupTogetherGrouping, 10000.}],
Cell[438574, 8787, 178, 3, 16, "Text",
 CellGroupingRules->{GroupTogetherGrouping, 10000.}],
Cell[438755, 8792, 11565, 332, 915, "Input",
 CellGroupingRules->{GroupTogetherGrouping, 10000.}],
Cell[450323, 9126, 841, 28, 13, "Output",
 CellGroupingRules->{GroupTogetherGrouping, 10000.}],
Cell[451167, 9156, 1841, 38, 62, "Output",
 CellGroupingRules->{GroupTogetherGrouping, 10000.}],
Cell[453011, 9196, 1136, 33, 31, "Output",
 CellGroupingRules->{GroupTogetherGrouping, 10000.}],
Cell[454150, 9231, 1075, 28, 31, "Output",
 CellGroupingRules->{GroupTogetherGrouping, 10000.}],
Cell[455228, 9261, 1075, 28, 31, "Output",
 CellGroupingRules->{GroupTogetherGrouping, 10000.}],
Cell[456306, 9291, 1504, 23, 187, "Output",
 CellGroupingRules->{GroupTogetherGrouping, 10000.}]
}, Open  ]],
Cell[457825, 9317, 1592, 49, 116, "Input"]
}, Open  ]]
}
]
*)

(* End of internal cache information *)

(* NotebookSignature GwTumdu9zVdaxCKT3aYXtlHx *)
