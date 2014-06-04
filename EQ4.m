(* ::Package:: *)

(*El paquete EQ4 continene las funciones necesarias para la construcci\[OAcute]n de elementos
de cuatro nodos con dos grados de libertad por cada nodo. Se emplementan elemenos 
isoparm\[EAcute]tricos para la construcci\[OAcute]n de los elementos. Este paquete tambi\[EAcute]n incluye
funciones para la generaci\[OAcute]n de las mallas a utilizar con los elementos Q4.*)

EQ4[Coor_,Elem_,Y_,\[Nu]_]:=Module[{N1,N2,N3,N4,xcoor,ycoor,x1,y1,x2,y2,x3,y3,x4,y4,x,y,i},(
xcoor=Table[0,{i,4},{i,1}];
ycoor=Table[0,{i,4},{i,1}];
For[i=1,i<5,i++,(
xcoor[[i]]=Coor[[Elem[[i]],1]];
ycoor[[i]]=Coor[[Elem[[i]],2]]
)];
x1=xcoor[[1]];x2=xcoor[[2]];x3=xcoor[[3]];x4=xcoor[[4]];
y1=ycoor[[1]];y2=ycoor[[2]];y3=ycoor[[3]];y4=ycoor[[4]];
N1=1/4(1-\[Xi])(1-\[Eta]);
N2=1/4(1+\[Xi])(1-\[Eta]);
N3=1/4(1+\[Xi])(1+\[Eta]);
N4=1/4(1-\[Xi])(1+\[Eta]);
x=Simplify[N1 x1+N2 x2+N3 x3+N4 x4];
y=Simplify[N1 y1+N2 y2+N3 y3+N4 y4];
JmatQ4={{D[x,\[Xi]],D[x,\[Eta]]},{D[y,\[Xi]],D[y,\[Eta]]}};
jacQ4=Det[JmatQ4];
JinvQ4=Inverse[JmatQ4];
d\[Xi]dx=JinvQ4[[1,1]];
d\[Xi]dy=JinvQ4[[1,2]];
d\[Eta]dx=JinvQ4[[2,1]];
d\[Eta]dy=JinvQ4[[2,2]];
dN1dx=D[N1,\[Xi]]d\[Xi]dx+D[N1,\[Eta]]d\[Eta]dx;
dN1dy=D[N1,\[Xi]]d\[Xi]dy+D[N1,\[Eta]]d\[Eta]dy;
dN2dx=D[N2,\[Xi]]d\[Xi]dx+D[N2,\[Eta]]d\[Eta]dx;
dN2dy=D[N2,\[Xi]]d\[Xi]dy+D[N2,\[Eta]]d\[Eta]dy;
dN3dx=D[N3,\[Xi]]d\[Xi]dx+D[N3,\[Eta]]d\[Eta]dx;
dN3dy=D[N3,\[Xi]]d\[Xi]dy+D[N3,\[Eta]]d\[Eta]dy;
dN4dx=D[N4,\[Xi]]d\[Xi]dx+D[N4,\[Eta]]d\[Eta]dx;
dN4dy=D[N4,\[Xi]]d\[Xi]dy+D[N4,\[Eta]]d\[Eta]dy;
BmatQ4={{dN1dx,0,dN2dx,0,dN3dx,0,dN4dx,0},
{0,dN1dy,0,dN2dy,0,dN3dy,0,dN4dy},
{dN1dy,dN1dx,dN2dy,dN2dx,dN3dy,dN3dx,dN4dy,dN4dx}};
Cmat=Y/(1-\[Nu]^2){{1,\[Nu],0},{\[Nu],1,0},{0,0,(1-\[Nu])/2}};
BDBe=(Transpose[BmatQ4].Cmat.BmatQ4)jacQ4;

(*REalizamos la integraci\[OAcute]n dependiendo del n\[UAcute]mero de puntos de Gauss 
Especificados*)

intg=BDBe;
q1=-0.5773502691;q2=-q1;w1=1.0;w2=1.0;
(ke=(intg w1 w1/.{\[Xi]->q1,\[Eta]->q1})+(intg w1 w2/.{\[Xi]->q1,\[Eta]->q2})+
(intg w2 w1/.{\[Xi]->q2,\[Eta]->q1})+(intg w2 w2/.{\[Xi]->q2,\[Eta]->q2}))//MatrixForm
)];

MallaQ4[XNodos_,YNodos_,bboxL_,bboxH_]:=Module[{res,i,k,kk},(
NNodos=XNodos YNodos;
NElementos=(XNodos-1)(YNodos-1);
p0=Table[0,{i,XNodos YNodos},{i,2}];
i=1;
Do[
Do[
p0[[i,1]]=(bboxL/(XNodos-1))(n-1);
p0[[i,2]]=(bboxH/(YNodos-1))(yn-1);i++,
{n,XNodos}]
,{yn,YNodos}];
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
,{yn,YNodos-1}];
Graphics[GraphicsComplex[p0,{Line[Elementos],Red,PointSize[Large],Point[p0]}]]
)]

(*FUNCION PARA REALIZAR MALLADOS NO PROPORCIONALES*)

MallaQ42[XNodos_,YNodos_,p0_]:=Module[{i,k,kk},(
NNodos=XNodos YNodos;
NElementos=(XNodos-1)(YNodos-1);
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
,{yn,YNodos-1}];
Graphics[GraphicsComplex[p0,{Line[Elementos],Blue,PointSize[Large],Point[p0]}]]
)]


(*FUNCI\[CapitalOAcute]N PARA GRAFICAR LA DEFORMADA*)

MallaDeform[XNodos_,YNodos_,p0_,p1d_]:=Module[{i,k,kk},(
NNodos=XNodos YNodos;
NElementos=(XNodos-1)(YNodos-1);
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
,{yn,YNodos-1}];
Deform=Graphics[{GraphicsComplex[p0,{Line[Elementos],Blue,PointSize[Large],Point[p0]}],GraphicsComplex[p1d,{Line[Elementos],Red,PointSize[Large],Point[p1d]}]}]
)]

(*FUNCION PARA CALCULAR LA MATRIZ DE RIGIDEZ DE TODOS LOS ELEMENTOS:*)
SMQ4[p_]:=Module[{i},(
MBDBe=Array[0&,{NElementos,8,8}];
For[i=1,(i<NElementos+1),i++,(
EQ4[p,Elementos[[i,All]],10*10^6,0.33];
MBDBe[[i]]=ke;
)]
)]

GlobalMK[p0_]:=Module[{i,j,k,l,posrow,pos,poscel,posfel},(
SMQ4[p0];
GlobalK=Array[0&,{NNodos 2,2 NNodos}];

For[i=1,i<(NNodos+1),i++,(
	For[j=1,j<NElementos+1,j++,(
		For[k=1,k<5,k++,(
			If[Elementos[[j]][[k]]==i,(
				For[l=1,l<5,l++,(
					pos=(Elementos[[j]][[l]])2;
					posrow=2 i;
					poscel=2*l;
					posfel=2*k;
					GlobalK[[posrow-1,pos-1]]=MBDBe[[j]][[posfel-1,poscel-1]]+GlobalK[[posrow-1,pos-1]];
					GlobalK[[posrow-1,pos]]=MBDBe[[j]][[posfel-1,poscel]]+GlobalK[[posrow-1,pos]];
					GlobalK[[posrow,pos-1]]=MBDBe[[j]][[posfel,poscel-1]]+GlobalK[[posrow,pos-1]];
					GlobalK[[posrow,pos]]=MBDBe[[j]][[posfel,poscel]]+GlobalK[[posrow,pos]];
				)]
			)]
		)]
	)]
)]
)]
