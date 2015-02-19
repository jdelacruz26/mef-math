(* ::Package:: *)

MallaQ4[XNodos_,YNodos_,bboxL_,bboxH_]:=(
XNodos=3;YNodos=2;bboxL=10.0;bboxH=5.0;
NNodos=XNodos* YNodos;
NElementos=(XNodos-1)(YNodos-1);
p0=Table[0,{i,XNodos*YNodos},{i,2}];
i=1;
Do[
Do[
p0[[i,1]]=(bboxL/(XNodos-1))(n-1);
p0[[i,2]]=(bboxH/(YNodos-1))(yn-1);i++,
{n,XNodos}]
,{yn,YNodos}]
Elementos=Table[0,{i,NElementos},{i,5}];
k=1;
kk=1;
Do[
For[i=1,(i<(XNodos+1)),i++,(
res=Round[yn  XNodos]-k;
If[res>0,
Elementos[[kk,All]]=Round[{k,k+1,k+(XNodos+1),k+XNodos,k}];
kk++
];
k++
)];
,{yn,YNodos-1}]
Graphics[GraphicsComplex[p0,{Line[Elementos],Red,PointSize[Large],Point[p0]}]]
)
