(* ::Package:: *)

(* ::Title:: *)
(*Droplet Optics*)


(* ::Chapter:: *)
(*functions*)


(* ::Section:: *)
(*compiled base functions*)


(* ::Subsection::Closed:: *)
(*directorAnalyticalC*)


(* ::Subsubsection:: *)
(*compiled code*)


{\[Theta]Center,\[Theta]Radius}={\[Pi]/2,(*0*)\[Pi]/2};
With[{\[Theta]Center=\[Theta]Center,\[Theta]Radius=\[Theta]Radius,\[Alpha]=\[Alpha],rd=rd},
directorAnalyticalC=Compile[{{x,_Real,0},{y,_Real,0},{z,_Real,0}},
	Module[{r,\[Theta]Bottom,\[Theta]Top,h,hr,\[Beta],\[Theta]zExp},
		r=Sqrt[x*x+y*y];
		\[Theta]Bottom=\[Theta]Center-((\[Theta]Center-\[Theta]Radius) Exp[-\[Alpha] (rd-r)]+\[Theta]Radius);
		h=Sqrt[rd*rd-r*r];
		hr=-r/h;
		\[Beta]=10.0/h;
		\[Theta]Top=0.5\[Pi]-ArcTan[-hr];
		\[Theta]zExp=(\[Theta]Top-\[Theta]Bottom)Exp[-\[Beta](h-z)]+\[Theta]Bottom;
		erC[\[Theta]zExp,If[x==0.0\[And]y==0.0,0.0,\[Pi]+ArcTan[x,y]]]
	]
	,CompilationTarget->"C",RuntimeOptions->"Speed",Parallelization->True,CompilationOptions->{"InlineCompiledFunctions"->True,"InlineExternalDefinitions" -> True}
];
]


(* ::Subsection::Closed:: *)
(*rotaMatrix\[Theta]\[Phi]C, rotaMatrix\[Theta]\[Phi]InvC*)


(* ::Text:: *)
(*same result as ...*)


(* ::Input:: *)
(*rotaMatrix=rotaMatrix\[Theta]\[Phi]C[\[Theta],\[Phi]]*)
(*rotaMatrixInv=rotaMatrix\[Theta]\[Phi]InvC[\[Theta],\[Phi]]*)


(* ::Subsubsection::Closed:: *)
(*derivation*)


(* ::Input:: *)
(*RotationMatrix[\[Pi]/2+\[Phi],{0,0,1}] . RotationMatrix[\[Theta],{1,0,0}] . {0,0,1};*)
(*RotationMatrix[\[Pi]/2+\[Phi],{0,0,1}] . RotationMatrix[\[Theta],{1,0,0}]*)
(*RotationMatrix[-\[Theta],{1,0,0}] . RotationMatrix[-(\[Pi]/2+\[Phi]),{0,0,1}]*)
(*(** test **)*)
(*RotationMatrix[\[Pi]/2+\[Phi],{0,0,1}] . RotationMatrix[\[Theta],{1,0,0}] .*)
(*RotationMatrix[-\[Theta],{1,0,0}] . RotationMatrix[-(\[Pi]/2+\[Phi]),{0,0,1}]//Simplify//MatrixForm*)


(* ::Subsubsection:: *)
(*code*)


rotaMatrix\[Theta]\[Phi]C=Compile[{{\[Theta],_Real,0},{\[Phi],_Real,0}},
	(** hardcoded output of: RotationMatrix[\[Pi]/2+\[Phi],{0,0,1}].RotationMatrix[\[Theta],{1,0,0}] **)
	{{-Sin[\[Phi]],-Cos[\[Theta]] Cos[\[Phi]],Cos[\[Phi]] Sin[\[Theta]]},{Cos[\[Phi]],-Cos[\[Theta]] Sin[\[Phi]],Sin[\[Theta]] Sin[\[Phi]]},{0,Sin[\[Theta]],Cos[\[Theta]]}}
	,CompilationTarget->"C",RuntimeOptions->"Speed",Parallelization->True
];

rotaMatrix\[Theta]\[Phi]InvC=Compile[{{\[Theta],_Real,0},{\[Phi],_Real,0}},
	(** hardcoded output of: RotationMatrix[-\[Theta],{1,0,0}].RotationMatrix[-(\[Pi]/2+\[Phi]),{0,0,1}] **)
	{{-Sin[\[Phi]],Cos[\[Phi]],0},{-Cos[\[Theta]] Cos[\[Phi]],-Cos[\[Theta]] Sin[\[Phi]],Sin[\[Theta]]},{Cos[\[Phi]] Sin[\[Theta]],Sin[\[Theta]] Sin[\[Phi]],Cos[\[Theta]]}}
	,CompilationTarget->"C",RuntimeOptions->"Speed",Parallelization->True
];


(* ::Subsection::Closed:: *)
(*spherical base vectors, er, erC*)


er[\[Theta]_,\[Phi]_]:={Sin[\[Theta]]Cos[\[Phi]],Sin[\[Theta]]Sin[\[Phi]],Cos[\[Theta]]};
e\[Theta][\[Theta]_,\[Phi]_]:={Cos[\[Theta]]Cos[\[Phi]],Cos[\[Theta]]Sin[\[Phi]],-Sin[\[Theta]]};
e\[Phi][\[Theta]_,\[Phi]_]:={-Sin[\[Phi]],Cos[\[Phi]],0};


(* ::Subsubsection:: *)
(*compiled code*)


erC=Compile[{{\[Theta],_Real,0},{\[Phi],_Real,0}},
	{Sin[\[Theta]]Cos[\[Phi]],Sin[\[Theta]]Sin[\[Phi]],Cos[\[Theta]]}
	,CompilationTarget->"C",RuntimeOptions->"Speed",Parallelization->True
];


(* ::Section:: *)
(*high-level functions*)


(* ::Subsection::Closed:: *)
(*director*)


(* ::Input:: *)
(*director2[{x,y,z}]*)


(* ::Subsubsection::Closed:: *)
(*debug: timing*)


(* ::Input:: *)
(*AbsoluteTiming[Do[Normalize[{u[x,y,z],v[x,y,z],w[x,y,z]}/.solution];,{i,10000}]][[1]]*)
(*AbsoluteTiming[Do[Normalize[{u[x,y,z]/.solution,v[x,y,z]/.solution,w[x,y,z]/.solution}];,{i,10000}]][[1]]*)
(*AbsoluteTiming[Do[Normalize[{u[x,y,z],v[x,y,z],w[x,y,z]}]/.solution;,{i,10000}]][[1]]*)
(*AbsoluteTiming[Do[director[{x,y,z}];,{i,10000}]][[1]]*)


(* ::Subsubsection:: *)
(*code*)


(* ::Input:: *)
(*director2[{x_,y_,z_}]:=ComplexExpand[Normalize[{u[x,y,z],v[x,y,z],w[x,y,z]}/.solution]];*)


director[{x_?NumericQ,y_?NumericQ,z_?NumericQ}]:=ComplexExpand[Normalize[{u[x,y,z],v[x,y,z],w[x,y,z]}/.solution]];


(* ::Subsection::Closed:: *)
(*tiltDirector*)


(* ::Text:: *)
(*same as Subscript[e, r]*)


tiltDirector[\[Theta]_,\[Phi]_]:={Cos[\[Phi]] Sin[\[Theta]], Sin[\[Phi]]Sin[\[Theta]], Cos[\[Theta]]};


(* ::Subsection::Closed:: *)
(*calcOseenFrankEnergy*)


(* ::Text:: *)
(*source: Borthagaray, Walker, "The Q-tensor model with uniaxial constraint"  (Elsevier, 2021), Eq. 10*)


calcOseenFrankEnergy[solution_]:=Module[{k1,k2,k3,energy},
	{k1,k2,k3}={1.0,1.0,1.0}; (** Frank elastic LC constants for 5CB **)
	energy=0.5 NIntegrate[
		k1 Div[director[{x,y,z}],{x,y,z}]^2 
		+ k2 (director[{x,y,z}] . Curl[director[{x,y,z}],{x,y,z}])^2
		+ k3 Norm[Cross[director[{x,y,z}],Curl[director[{x,y,z}],{x,y,z}]]]^2
		,{x,y,z}\[Element]halfBall];
	
	energy
]


(* ::Input:: *)
(*Grad[]*)


(* ::Input:: *)
(*Image[CompressedData["*)
(*1:eJztnV1MU2kax8/u3uzFudheckHCpmnSCxJCSJMtIZP2ouyYJmA0pAENKWRc*)
(*IONMi66gxoJZKGYsGaVxhB1tHBsmW3ZoyFoNXYcq1A04Q52pAaJVwLTKGUFR*)
(*mkpaOfRlz1e/v9sjUHx/iQZKz2ff/3me93me9+mfGxUH/vZ7BEGO/5H470BD*)
(*u/jYsYZTB/9E/FLz5fGWpi+PfLbvy78faTpy7C+NfyBe/Cvxr+B3CEL+vAWB*)
(*QCAQCAQCgUAgEAgEAoFAIBAIBAKBQCAQCAQCgUAgkHxi463Doms7qXd4Wd0t*)
(*wLH/9TdVoAiCcA+qLc9xVvcOgewmaBFVcxGy0FLHqpSAx6at4gnaRh47RjvE*)
(*hQjnuHl1k8X9QyC7BrDxfOaB042//7WvnMO2lNz2vv0o/4xllfjhU1KpiEhj*)
(*87C3fwhk9+F/pKssYFdKYMnYwEGLuyfXt8C6faCay+HJBx1ewNb+IZDdCPtS*)
(*ejfXfwBBhCrrKks7hEDyAdalRHuMcHIE+dhgWUrgvV1bTsyNKnUOPxv7g0Dy*)
(*BZaltD7XX00usG6zuNnYHQSSN7ArJfBUX12IIGhpnx26d5CPC3altGpWcMik*)
(*bL3RxcLeIJA8gk0pgXfWzsK9kEXyYZNXFGIeeSloubxvHMNhJB+SCjalxEyU*)
(*EOSQftHHwrntDJvu6QuV4s8v6n8Y1vc2i4hnQ6Gk78H6Tp8WZLfDopTAM4Os*)
(*iFRSXkfC/U/0h1XmZaZmEKyOtfNRpLDT+g4aJkhSWJSS29JWQBklrsaWt9Wr*)
(*wDV8Sje3EXph1aoSIkiTEcvbS4JsD+zV4G3O9VdQSkJkBtfeeYR7bBoRUqKx*)
(*+fbOJUHYBXhcD6cnbmjk1Oy6UNz23diU7aHLk+2I2VwxfYFSSuIozHunZogM*)
(*73OFmp/hXAmyXaxNdTNGKZ/9uyjAhkNXzW81YRup3xuB3/t2GcOw397mXMWL*)
(*v3Lcu6Vr+xRlbq6opj5AjYhaIIOgtQZnntlMytYTZ14mPUxeifyowZH27Nrr*)
(*MLRT139QxEX31nijYJKzpIGr1j/N9oOl7nC9EQvslBzJzBhKh/J2y3LGh8Zt*)
(*Gm6CRBhYMrdWNxsXM/ukiMF//4amhrD2dHk8C/gdusp4Fwi8Tou6uoSlo2wj*)
(*tJSKZIZnOTwDfIv6Q3tQSsGYA8JvMi1lu5coKRFjZdnSXk7tVamfnLYlZFyv*)
(*EKLVOsdG5p9MQin5lkyn67S2rFzet1PdxFCpUFlfZ7FxDH635VQB5TqblvGY*)
(*P3UfNb3IM6PESCnHVD6OGZv2npRCMYec8rMxUtoC67ZeYUqLsz7ZXVzcYMzK*)
(*zYkvpU3PzLXWs7ezTM4y5fGxIz87PPY+aQJHDqzPjt56km9GCUopIaGYQ275*)
(*2VgpEYPFaWzgkwMpodHZWDK2cIS9tvWshn0cKW16Zq+3fD60EDrc+sKt23Pv*)
(*09w/45dymk2Zu5vx2JjTVRcR7tB+/ZNQuT147/Xlb/E9lFIimFlkzvnZeFIK*)
(*GSZR99TbOBut/6wRlmRpkrZipQRwbLRdKG27+m8jw7Dh8rHaDIJ4752GRhTh*)
(*1hoWWVEScBpqyeeUtM8eNPebq5bL/7Tnb3VWPCl5Zoe75IJAeKXe6MSxSZ1K*)
(*LuKiaNl+heamwxM1rlJLCWC3WgUFqODkaMaBox0iWOeQ6zMirpRChonTYFyK*)
(*Hp34skkRaZL8njljV315wEo2GbF32NR1Vb2YixSUVSk0pkcR059IKYGVu53i*)
(*wpiAxoH+uXfJThxfdViudXf06Yf1A4M/XG0qRtA6nYMSn985ekwSip3w6rq/*)
(*/wnzux2mi1RVEl+qHLR7ktgX6gKJDYt7pphr9HuXxtQHNVNxrDDAV341alrE*)
(*XL6o5vBhaRmnrPaM7u4iMwiBd/GOTiWj0h884g2BUGCdtKxge4uQE1il0EAq*)
(*Onj+Uqf8aJems0lEfxwov31sNeKKU0op6CzlMn/fXnxT3QElxRFCBiSQUsgw*)
(*faqxrUX8hXR+SmOf/8BlkDEnJDt/9ay8uVPT1SRC4wX6kkTw0gLgb3653tLY*)
(*aX5GBb7x5Zsn+EUoUq61Bx1CfN4g5zOnIxmYY/zGVauqUjZgTxXWeG1VUdNQ*)
(*JgxOD3uUr5qIWRFG+KXfyXkcnvzb6RUf88qCSSUuRMVdY1jA635nVRUi0fUb*)
(*4LW1Y9/OS2nLZaxnAvx8xQ/z5BMgsJg0jsOThoPnfWpSK5Wa0QUmKcFskgE5*)
(*jeeMAYv66sCRi7qncqhkTSSl4HwhyjARfs4ZPtnaKMalxIz1zBkJFUOPqeHK*)
(*TN6RqCRyblKirBhfornvZs6KjrZFhcHDo/plVHQdeGe/rW0hnJhUPiAZUSG2*)
(*LG6mw3TA+xabMSr3NcdE7QB2o4WHIoIu65vwuwE25gdrOQgq0dppK0Zeb4yU*)
(*yKr+ns92kZTC/hT6KKPOOYu5EmGX740YM2DEushug8jkpxdYW4HkllTaSiYl*)
(*WjVRhonUV0l1RMlcgDj3P+yJFH6IXKS08cTQUIpKLky7g6OXnijFhsGphmb0*)
(*0XkK4+P7/XWnTUspnzqBZzLaaHC+D7zodVzri4kN0saroFL3KMZZfGVpK0OC*)
(*c7f4UtpmdkpKbJKmWctkl0ymjDYBuTUaSiKlkGEqUJopn5ny+uKapK3tkdK6*)
(*Q1eHIvyIiAd4YWoujhsGB2/GOwQc6vBoEXffKQuWxjPHt6A/jEZHL72uB3Ov*)
(*ojZm3LZP++yxc7rAlKFaT2opjpSIN5zTxtmQOXF8ZXrwZGP9iX90n5AJhJ9d*)
(*ND8NOqXgzT1127WJm+oqsvCAX6WZfBM8MbA2N9RRXyOTlhVypaeH5tbCTnkv*)
(*SOkDQNdO0zQYXLmESpJKKWiY6Bk9GYsolfTPxD/eh5cSs/4iUsv0iwkqeYKG*)
(*lUB0Mh0p0cJEOOV9v76P+wbc56VcxEA5RHxbw6in4JTF7Y8jJfCbWXkwYTbQ*)
(*+7C/SiA3zFPv9mHmMwJE2GpeAviy3aiRCwhnVqLUT6+sY1Y14T9X98/Rjq1n*)
(*tv+QRPMTITrgnlQLOUhEaAhKKR7h4bsi9VROFdTJpUQYppl+CTXpbr/9YrpX*)
(*yGkxLiVQ7geXEhNYi0we0XYqMgEUzuZMfwWHOQd+uzlVAhcsm5rJtyeqmvA4*)
(*Br8epOKEjD6SS4n+K/NLKIJ3WFqGJkysU5cZPt+nV+VQugBrlvbC0M2kdszc*)
(*SaqtaHC5AV2fWdZmeRU88x2S0u4OOwQ+J9JcNJlWctpXKinR2VjySALRJ7xk*)
(*BdsfXEp0aVDE3AQsj7YSdioiARRxdQ59W+vVa11MsB3lNScPO4D1qR7CJiWo*)
(*mvBhY12SOsb85WCVNj3Ou32yqgRSouZZEYP1HdW8mpIJvas4UqJcSlTc1HVe*)
(*Q9LTVlNCXEZF/0xAkNslpegI3u4OO6yYmpjpdPi9yo6UUqITstSDPYlJ2toJ*)
(*KeHO/w70nJAWUQkg3DN7c+RheLgaX5k439h73w02PTNXZIxpKm0wPAldA1iL*)
(*d4iABEJv875yTA6r63jBsN5WjnMl4s50JZgrUWM7wtmgbyNlxRJKiZJbMhcl*)
(*ZV4pFL8Ky2tETR9S55VWzcepO12iML9McCa7iMDzEIm04NmRhpQYw8RJvoYo*)
(*7P4HC5k2XIYG5rXwxYnZOngr5jay+Ftybmz+hcsxoe/51vpsnAxcC1XDt6/3*)
(*Dj0Ke5pRGZ/KrwOTBdoPpB8/Dfp5+jpIGxQ8fRyzjVxWBhJhwcUIoYUVJBHT*)
(*tA8UwaMnwuFjmB7A1F1NLqVk60kTVDsw6WMKXt25iRc4drdXHsy2F4qaL1ux*)
(*4KwxpZQA7vyPkqx2aLvh3MaQdraE+Xe5d0dJR0qUYSqX6xcSHcvvmRtWkR4F*)
(*Mwx5st4JzINN9JFzZOY1cbN2AvMHLyCrCB6+NDmgFIeXtWw8HW75hCdWDkwu*)
(*hT5c/wvrpVYpF0V4MtW1Sczv9zhuXWwWB+sfUNHpURd5LWBlPONzCAHwJZOS*)
(*cC+FX025I/JKuHNIzkHQqm/sdM1DZlIiHuwnCsLN35bbppGg+/UL/q0kDh5l*)
(*DjhC9WQg3UY8HH7sOXsrEGqBNXixhHlNucYcttKVErvkWu2we0iv2oEpTUkQ*)
(*awUuU0sZwgtb7UjnzgJKJEP6wgOa6dWQgYsXdthaf9Anocqi2r4xGI3D18/J*)
(*JceMIdOwt6QE3It3rnd/sb8MRajn6oh9JQub4nXoahhXhYV16FBKOZKyBm9C*)
(*311HtyAQKbRDsXNq/zNjY2lU0RrwzAyp6kSSI6qvTjfJVYO2l+TQJfY1oBCQ*)
(*nphcM/KLc/7ugOIT0geoJ0cSIM3QeB/jmxUI5BfNC+6YvBJH2DG2nHWPQfDG*)
(*pqnaeSmBl5O9h3h0WVfQ/eYpjBm7lMGa8LheeqZ85FLaJc0DfYv64/GL8Fkj*)
(*uKicJrsF6Vlsyzr4svnMfpVpgXlY+b2ucW19KflIabmRTg4+jGD0sqJ7ai31*)
(*21PwcUsJ7IbuMoQ1uaPV3oENadNiY0539NJ05OSUqbfMdOEn00aPrTaSPpfl*)
(*isawvetw/E7LJa3B/kEfwnkE/vrJ0xWoo/QArtFvzM+jb1Z2DfGC8Tv4hUoQ*)
(*CA1tX0Lry9KDyYSiJZrpXeLoQyA7C1XuFbtQMQWB7ij5kU2GQD48VNEyr8Ww*)
(*kFHvmuBKpbz+6goIhD3IpUAVEfVg6W3GlOJU9M9FhxzIBMeI9mS9KFAJwhXJ*)
(*lJp/WRey7qQMgex61hcMimr1vTcZD3ImqRTzzbM+bOKCjMfhStt1d+Yp7VCp*)
(*Q7ICs1AcCsJDIHuJTfe0trEzu96JdFKJJ49onBVsJnAhMt5O4HUaFTyEI+gY*)
(*z1y2EMhuZtMzO3iy58cs83FMJDwqfk4vDSiO7eNBQofcIzoVQCD5DlmC3hmr*)
(*I3wVe51eDIGOhEdnlF6aFSWJ2x3Thiy4bBkCyXcS1IeAtdnvr99ZSafgAfim*)
(*1EVx2g7QHQwS9CalWlFl2ScfAtl1EDq63SEuDevPWR9Y5i8oTje1RNWEo4fj*)
(*rBsCv0107kOjF86QmzgNTRxek96Rvx16IZAQZMOlOL18aWLan4Zth2PjWqVS*)
(*bXpKTo3AoqGWlyilCzwPdfLysEaL1OZOYzNf2jH2fE8s1IJAsoaeBNGtBjbJ*)
(*8teC4GLqcPBXj568Ils23ddU1jBrxMiXF43NUqoxKfMGCORjhQ7NobyG7+2P*)
(*b2tqxAm+Ds9j0yi09reAavlYTq9c3vK7rWdLqB7FwGPTyi7u9LpHCGQnASv3*)
(*emWkYUIFjdHfBxGCztvypQp1/5WzNRy6hwDxokTwxVeX1cdqygp2et0jBJIX*)
(*hH3XEhJsxxH5IpQSBAKBQCAQCAQCgUAgEAgEAslD/g9HseW6*)
(*"], "Byte", ColorSpace -> "RGB", ImageSize -> {115.6, Automatic}, Interleaving -> True, MetaInformation -> Association["Comments" -> Association["Software" -> "Greenshot"]]]*)


(* ::Input:: *)
(*Subscript[sum, ij](\!\( *)
(*\*SubscriptBox[\(\[PartialD]\), \(i\)]*)
(*\*SubscriptBox[\(n\), \(j\)]\))(\!\( *)
(*\*SubscriptBox[\(\[PartialD]\), \(i\)]*)
(*\*SubscriptBox[\(n\), \(j\)]\))*)


(* ::Input:: *)
(*calcOseenFrankEnergyOneConstant[solution]*)


(* ::Input:: *)
(*u[x,y,z]/.solution*)
(*D[u[x,y,z],x]/.solution*)


calcOseenFrankEnergyOneConstant[solution_]:=Module[{k,energy},
	k=1.0; (** Frank elastic LC constants for 5CB **)
	energy=0.5 NIntegrate[
		((D[u[x,y,z],x]^2+D[u[x,y,z],y]^2+D[u[x,y,z],z]^2)
		+(D[v[x,y,z],x]^2+D[v[x,y,z],y]^2+D[v[x,y,z],z]^2)
		+(D[w[x,y,z],x]^2+D[w[x,y,z],y]^2+D[w[x,y,z],z]^2))/.solution
		,{x,y,z}\[Element]halfBall];
	
	energy
]


(* ::Input:: *)
(*solution*)


(* ::Input:: *)
(*divDirector=Div[ComplexExpand[Normalize[{u[x,y,z],v[x,y,z],w[x,y,z]}/.solution]],{x,y,z}]/.{x->0.1,y->0.1,z->0.1}*)


(* ::Input:: *)
(*k2Director=ComplexExpand[Normalize[{u[x,y,z],v[x,y,z],w[x,y,z]}/.solution]] . Curl[ComplexExpand[Normalize[{u[x,y,z],v[x,y,z],w[x,y,z]}/.solution]],{x,y,z}]/.{x->0.1,y->0.1,z->0.1}*)


(* ::Input:: *)
(*k3Director=ComplexExpand[Normalize[{u[x,y,z],v[x,y,z],w[x,y,z]}/.solution]] . Curl[ComplexExpand[Normalize[{u[x,y,z],v[x,y,z],w[x,y,z]}/.solution]],{x,y,z}]/.{x->0.1,y->0.1,z->0.1}*)


(* ::Input:: *)
(*divDirector=Div[ComplexExpand[Normalize[{u[x,y,z],v[x,y,z],w[x,y,z]}/.solution]],{x,y,z}]/.{x->0.1,y->0.1,z->0.1}*)


(* ::Input:: *)
(*Div[Normalize[{u[x,y,z],v[x,y,z],w[x,y,z]}],{x,y,z}]/.solution/.{x->0.1,y->0.1,z->0.1}*)


(* ::Input:: *)
(*Clear[f,g,h]*)
(*Div[{f[x,y,z],g[x,y,z],h[x,y,z]},{x,y,z}]*)


(* ::Text:: *)
(*source: Borthagaray, Walker, "The Q-tensor model with uniaxial constraint"  (Elsevier, 2021), Eq. 11*)


(* ::Input:: *)
(*calcOseenFrankEnergyOneConstant[solution_]:=Module[{k1,k2,k3,energy},*)
(*	{k1,k2,k3}={1.0,1.0,1.0}; (** Frank elastic LC constants for 5CB **)*)
(*	energy=0.5 NIntegrate[*)
(*		k1 Div[director[{x,y,z}],{x,y,z}]^2 *)
(*		+ k2 (director[{x,y,z}] . Curl[director[{x,y,z}],{x,y,z}])^2*)
(*		+ k3 Norm[Cross[director[{x,y,z}],Curl[director[{x,y,z}],{x,y,z}]]]^2*)
(*		,{x,y,z}\[Element]halfBall];*)
(*	*)
(*	energy*)
(*]*)


(* ::Input:: *)
(*solution=solveHemisphereLC[0.0,0.0];*)


(* ::Input:: *)
(*calcOseenFrankEnergy[solution]*)


(* ::Subsection::Closed:: *)
(*calcXYProjectedDirectorLength*)


calcXYProjectedDirectorLength[lcDirector_,waveVector_,nLCe_]:=Norm[lcDirector-waveVector(lcDirector . waveVector)]nLCe


(* ::Subsection::Closed:: *)
(*calcXYnLC*)


(* ::Text:: *)
(*calculate refractive index that a wave with a polarization vector going through liquid crystal with optical indicatrix feels*)
(*Note: This is NOT the extraordinary refractive index nLCe' needed for Jones matrix formalism*)


(* ::Input:: *)
(*{nLCe,nLCo}*)
(*lcDirector=Normalize@{1,0,1};*)
(*nLCMod=calcXYnLC[lcDirector,polarizationVector]*)


calcXYnLC[lcDirector_,polarizationVector_,nLCe_,nLCo_]:=Module[{rota,rotatedP},
rota=RotationMatrix[{lcDirector,{1.0,0,0}}];  (** transformation matrix, so that LC director aligns with x-axis for simpler calculation **)
rotatedP=rota . polarizationVector;
(** intersection of ellipsoid E with ray r from plugging in ray through origin into ellipsoid; solution of: ( (px/ne)^2 + (py/no)^2 + (pz/no)^2 ) t^2 = 1 **)
1.0/Sqrt[((rotatedP[[1]]/nLCe)^2+(rotatedP[[2]]/nLCo)^2+(rotatedP[[3]]/nLCo)^2)]
];


(* ::Subsection::Closed:: *)
(*calcXYnLCe*)


(* ::Text:: *)
(*intersection of 3d ellipsoid on xy-plane = ellipse*)
(*This is the extraordinary refractive index nLCe' needed for Jones matrix formalism*)
(**)
(*using input (Subscript[n, e],Subscript[n, o]) = (\[CapitalDelta]n+Subscript[n, o],Subscript[n, o])*)


(* ::Input:: *)
(*{nLCe,nLCo}*)
(*(*nLCeMod=calcXYnLCe[lcDirector,polarizationVector]*)*)


(* ::Subsubsection::Closed:: *)
(*derivation*)


(* ::Input:: *)
(*(** project ellipsoid on plane **)*)
(*(** boundary points are identified as those whose tangent plane contains z-axis < = > normal vector has no z-component **)*)
(*Clear[a,b,c,d,e,f]*)
(*ellipsoid[x_,y_,z_]:=a x^2 + b y^2+ c z^2+2(d x y + e x z+f y z)-1*)
(*ellipsoidCut[x_,y_]:=ellipsoid[x,y,0]//Simplify   (** cross section through plane **)*)
(*ellipsoidCut[x,y]*)


(* ::Input:: *)
(*{a,b,c,d,e,f}=N@{3,1,5,1,2,2};*)
(*Grid[{*)
(*Show[*)
(*DiscretizeRegion[ImplicitRegion[ellipsoid[x,y,z]==0,{x,y,z}],Axes->True,Boxed->True,AxesLabel->{"x","y","z"},AxesStyle->Directive[Black,18],BoxStyle->Directive[Black]],*)
(*(*DiscretizeRegion[ImplicitRegion[projectedEllipse[x,y]<0\[And]z==-1,{x,y,z}]],*)*)
(*PlotRange->{All,All,{-1,1}}*)
(*],*)
(*DiscretizeRegion[ImplicitRegion[ellipsoidCut[x,y]<0,{x,y}],Frame->True]*)
(*}]*)


(* ::Input:: *)
(*(** set \[Phi] to zero, so that ellipse is aligned with coordinate axes and nLce can be read off easily **)*)
(*Clear[nLCe,nLCo]*)
(*ellipsoidLC=(RotationMatrix[{{0,0,1},{Sin[\[Theta]]Cos[\[Phi]],Sin[\[Theta]]Sin[\[Phi]],Cos[\[Theta]]}}] . {x,y,z})\[Transpose] . {{nLCo,0,0},{0,nLCo,0},{0,0,nLCe}} . (RotationMatrix[{{0,0,1},{Sin[\[Theta]]Cos[\[Phi]],Sin[\[Theta]]Sin[\[Phi]],Cos[\[Theta]]}}] . {x,y,z})//ComplexExpand//FullSimplify;*)
(*ellipsoidCoeffs=Simplify[Coefficient[ellipsoidLC,#]]&/@{x^2,y^2,z^2,x y,x z,y z};*)
(*projectedEllipseCoeffs={*)
(*ellipsoidCoeffs[[1]],*)
(*ellipsoidCoeffs[[2]],*)
(*ellipsoidCoeffs[[4]]*)
(*}/.{\[Phi]->0}//Simplify*)


(* ::Input:: *)
(*f[z_]:=m z +b*)
(*Integrate[Cos[2f[z]],z]*)


(* ::Subsubsection:: *)
(*code*)


calcXYnLCe[\[Theta]_,nLCe_,nLCo_]:=1/2 (nLCe+nLCo-(nLCe-nLCo) Cos[2\[Theta]])


(* ::Subsubsection:: *)
(*compiled code*)


calcXYnLCeC=Compile[{{\[Theta],_Real,0},{nLCe,_Real,0},{nLCo,_Real,0}},
	0.5(nLCe+nLCo-(nLCe-nLCo) Cos[2\[Theta]])
	,CompilationTarget->"C",RuntimeOptions->"Speed",Parallelization->True
];


(* ::Subsection::Closed:: *)
(*calcXYnLCeShadow*)


(* ::Text:: *)
(*projection of 3d ellipsoid on plane = ellipse*)
(*This is NOT the extraordinary refractive index nLCe' needed for Jones matrix formalism*)


(* ::Subsubsection::Closed:: *)
(*derivation*)


(* ::Input:: *)
(*(** project ellipsoid on plane **)*)
(*(** boundary points are identified as those whose tangent plane contains z-axis < = > normal vector has no z-component **)*)
(*Clear[a,b,c,d,e,f]*)
(*ellipsoid[x_,y_,z_]:=a x^2 + b y^2+ c z^2+2(d x y + e x z+f y z)-1*)
(*normal[x_,y_,z_]:=Grad[ellipsoid[x,y,z],{x,y,z}]*)
(*zValsProjectionBoundary=Solve[normal[x,y,z][[3]]==0,z][[1]];  *)
(*projectedEllipse[x_,y_]:=ellipsoid[x,y,z]/.zValsProjectionBoundary//Simplify*)
(*projectedEllipse[x,y]//Collect[#,{x,y}]&   (** projection on plane **)*)
(*ellipsoid[x,y,0]  (** cross section **)*)
(**)


(* ::Input:: *)
(*{a,b,c,d,e,f}=N@{3,1,5,1,2,2};*)
(*Grid[{*)
(*Show[*)
(*DiscretizeRegion[ImplicitRegion[ellipsoid[x,y,z]==0,{x,y,z}],Axes->True,Boxed->True,AxesLabel->{"x","y","z"},AxesStyle->Directive[Black,18],BoxStyle->Directive[Black]],*)
(*(*DiscretizeRegion[ImplicitRegion[projectedEllipse[x,y]<0\[And]z==-1,{x,y,z}]],*)*)
(*PlotRange->{All,All,{-1,1}}*)
(*],*)
(*DiscretizeRegion[ImplicitRegion[projectedEllipse[x,y]<0,{x,y}],Frame->True]*)
(*}]*)


(* ::Input:: *)
(*(** set \[Phi] to zero, so that ellipse is aligned with coordinate axes and nLce can be read off easily **)*)
(*ellipsoidLC=(RotationMatrix[{{0,0,1},{Sin[\[Theta]]Cos[\[Phi]],Sin[\[Theta]]Sin[\[Phi]],Cos[\[Theta]]}}] . {x,y,z})\[Transpose] . {{nLCo,0,0},{0,nLCo,0},{0,0,nLCe}} . (RotationMatrix[{{0,0,1},{Sin[\[Theta]]Cos[\[Phi]],Sin[\[Theta]]Sin[\[Phi]],Cos[\[Theta]]}}] . {x,y,z})//ComplexExpand//FullSimplify;*)
(*ellipsoidCoeffs=Simplify[Coefficient[ellipsoidLC,#]]&/@{x^2,y^2,z^2,x y,x z,y z};*)
(*projectedEllipseCoeffs={ellipsoidCoeffs[[1]]-ellipsoidCoeffs[[5]]^2/ellipsoidCoeffs[[3]],ellipsoidCoeffs[[2]]-ellipsoidCoeffs[[6]]^2/ellipsoidCoeffs[[3]],ellipsoidCoeffs[[4]]-ellipsoidCoeffs[[5]]ellipsoidCoeffs[[6]]/ellipsoidCoeffs[[3]]}/.{(*nLCe->1.7,nLCo->1.5,\[Theta]->0.1,*)\[Phi]->0}//Simplify*)


(* ::Subsubsection:: *)
(*code*)


calcXYnLCeShadow[\[Theta]_,nLCe_,nLCo_]:=1/4 (nLCe+3 nLCo+(-nLCe+nLCo) Cos[2 \[Theta]]+2 (nLCe-nLCo) Sin[\[Theta]]^2-(8 (nLCe-nLCo)^2 Sin[2\[Theta]]^2)/(nLCe+nLCo+(nLCe-nLCo) Cos[2\[Theta]]))


(* ::Subsection::Closed:: *)
(*calcLCrefIndicesXY*)


(* ::Input:: *)
(*{nLCe,nLCo}*)
(*lcDirector={0,0,-1};  (** both n's are nLCo **)*)
(*lcDirector=Normalize[{1,0,0}];  (** nLCx is nLCe, nLCy is nLCo **)*)
(*lcDirector=Normalize[{0,1,0}];  (** nLCx is nLCo, nLCy is nLCe **)*)
(*lcDirector=Normalize[{1,1,0}];  (** nLCx = nLCy is mix of {nLCe,nLCo} (larger than mean, since it lies on ellipsoid) **)*)
(*lcDirector=Normalize[{1,1,1}];  (** nLCx = nLCy is smaller mix of {nLCe,nLCo} **)*)
(*xPolarizationVector={1,0,0};*)
(*yPolarizationVector={0,1,0};*)
(*{nLCx,nLCy}=calcLCrefIndicesXY[lcDirector,xPolarizationVector,yPolarizationVector]*)


calcLCrefIndicesXY[lcDirector_,xPolarizationVector_,yPolarizationVector_,{nLCe_,nLCo_}]:=Module[{rota,xRotated,yRotated,nLCx,nLCy},
	rota=RotationMatrix[{lcDirector,{1.0,0.0,0.0}}];  (** transformation matrix, so that LC director aligns with x-axis for simpler calculation **)
	xRotated=rota . xPolarizationVector;
	yRotated=rota . yPolarizationVector;
	(** intersection of ellipsoid E with ray R from plugging in ray through origin into ellipsoid: ( (px/ne)^2 + (py/no)^2 + (pz/no)^2 ) t^2 = 1 **)
	nLCx=((xRotated[[1]]/nLCe)^2+(xRotated[[2]]/nLCo)^2+(xRotated[[3]]/nLCo)^2)^(-0.5);
	nLCy=((yRotated[[1]]/nLCe)^2+(yRotated[[2]]/nLCo)^2+(yRotated[[3]]/nLCo)^2)^(-0.5);
	{nLCx,nLCy}
];


(* ::Subsection:: *)
(*calcLCn\[Phi]*)


(* ::Input:: *)
(*data=calcLCn\[Phi][tiltAngleTheta,tiltAnglePhi,xVals,yVals,zVals,compiledQ]*)


(* ::Text:: *)
(** LC director projected on xy plane angle = marker for liquid type (LC,FC,water)*)
(** compiled code is 20x faster than uncompiled code and gives the same results*)


(* ::Subsubsection::Closed:: *)
(*debug*)


(* ::Input:: *)
(*solution=solveHemisphereLC[0.0,0.0];*)
(*setDefaultParameters[];*)
(*\[Lambda]=4.7`*^-7;*)
(*setRefractionIndices[\[Lambda]];*)
(*{xVals,yVals,zVals,nxJones,nyJones,nzJones,d}=setJonesCalculusDiscretization[rDropletDim];*)
(*compiledQ=True;*)


(* ::Input:: *)
(*{tiltAngleTheta,tiltAnglePhi}={0.0,0.0};*)
(*directorAnglesRefIndices3d=calcLCn\[Phi][tiltAngleTheta,tiltAnglePhi,xVals,yVals,zVals,compiledQ]*)


(* ::Subsubsection:: *)
(*code*)


calcLCn\[Phi][tiltAngleTheta_,tiltAnglePhi_,xVals_,yVals_,zVals_,compiledQ_:True]:=Module[{
director,rotaMatrix,rotaMatrixInv,r,position,rotatedPosition,rotatedDirector,lcAngleTheta,lcAnglePhi,nLCeMod},

If[compiledQ,
(** compiled version **)
calcLCn\[Phi]C[tiltAngleTheta,tiltAnglePhi,xVals,yVals,zVals]
,
(** uncompiled version **)
director=tiltDirector[tiltAngleTheta,tiltAnglePhi];
rotaMatrix=RotationMatrix[{{0,0,1},director}];
rotaMatrixInv=RotationMatrix[{director,{0,0,1}}];

Table[
	position={x,y,z};
	r=Norm[position];
	(*{directorAngle,refIndexE}=*)If[r<=0.99,
		(** inside droplet **)
		If[position . director>0.0,
			rotatedPosition=rotaMatrixInv . position;
			rotatedDirector=rotaMatrix . director[rotatedPosition];
			(*rotatedDirector=rotaMatrix.directorAnalytical[rotatedPosition];*)
			lcAngleTheta=(*\[Pi]/2.0*)ArcCos[rotatedDirector[[3]]];(** how much is the director tilted away from Subscript[e, z]? (\[Theta]) **)
			lcAnglePhi=ArcTan[rotatedDirector[[1]],rotatedDirector[[2]]];   (** how much is the director turned in \[Phi]? **)
			nLCeMod=calcXYnLCe[lcAngleTheta,\[CapitalDelta]nLC+nLCo,nLCo];
			(*nLCeMod=calcXYnLCeShadow[lcAngleTheta,nLCe,nLCo];*)
			{lcAnglePhi,nLCeMod}
			,
			{noAngle,nOil}
		]
		,
		(** outside droplet **)
		{noAngle,nWater}
	]
,{x,xVals},{y,yVals},{z,zVals}]
]
]


(* ::Subsubsection::Closed:: *)
(*compiled code*)


Clear[x,y,z];
With[{\[CapitalDelta]nLC=\[CapitalDelta]nLC,nLCo=nLCo,noAngle=noAngle,nOil=nOil,nWater=nWater},
calcLCn\[Phi]C=Compile[{{tiltAngleTheta,_Real,0},{tiltAnglePhi,_Real,0},{xVals,_Real,1},{yVals,_Real,1},{zVals,_Real,1}},
Module[{r2,rotatedPosition,rotatedDirector,lcAngleTheta,lcAnglePhi,nLCeMod,rotaMatrix,rotaMatrixInv},
rotaMatrix=rotaMatrix\[Theta]\[Phi]C[tiltAngleTheta,tiltAnglePhi];
rotaMatrixInv=rotaMatrix\[Theta]\[Phi]InvC[tiltAngleTheta,tiltAnglePhi];
(*Parallel*)Table[
r2=x*x+y*y+z*z;
If[r2<=0.99*0.99,
If[{x,y,z} . erC[tiltAngleTheta,tiltAnglePhi]>0.0,   (** vector e_r is normal to interface plane **)
rotatedPosition=rotaMatrixInv . {x,y,z};
rotatedDirector=rotaMatrix . directorAnalyticalC[Compile`GetElement[rotatedPosition,1],Compile`GetElement[rotatedPosition,2],Compile`GetElement[rotatedPosition,3]];(*rotaMatrix.directorAnalyticalC[rotatedPosition\[LeftDoubleBracket]1\[RightDoubleBracket],rotatedPosition\[LeftDoubleBracket]2\[RightDoubleBracket],rotatedPosition\[LeftDoubleBracket]3\[RightDoubleBracket]];*)
lcAngleTheta=ArcCos[Compile`GetElement[rotatedDirector,3]];
lcAnglePhi=ArcTan[Compile`GetElement[rotatedDirector,1],Compile`GetElement[rotatedDirector,2]];   (** how much is the director turned in \[Phi]? **)
nLCeMod=calcXYnLCeC[lcAngleTheta,\[CapitalDelta]nLC+nLCo,nLCo];
{lcAnglePhi,nLCeMod},{noAngle,nOil}
]
,{noAngle,nWater}]
,{x,xVals},{y,yVals},{z,zVals}]
]
,CompilationTarget->"C",RuntimeOptions->"Speed",Parallelization->True,CompilationOptions->{"InlineCompiledFunctions"->True,"InlineExternalDefinitions"->True}]
];


(* ::Subsection::Closed:: *)
(*jinc*)


(* ::Text:: *)
(*jinc for microscope blurring kernel*)


jinc[x_,y_,r0_,x0_,y0_]:=Module[{rxy},
	rxy=Sqrt[(x-x0)^2+(y-y0)^2];
	r0 BesselJ[1, 2\[Pi] r0 rxy]/rxy
];


(* ::Subsection::Closed:: *)
(*dispersion relations n(\[Lambda])*)


(* ::Text:: *)
(*sources: *)
(*https://www.guidechem.com/trade/hfe-7200-ethyl-nonafluorobutyl-id6163721.html*)
(*https://refractiveindex.info/?shelf=other&book=5CB&page=Tkachenko-o*)
(*https://refractiveindex.info/?shelf=other&book=5CB&page=Tkachenko-e*)


(* ::Subsubsection:: *)
(*code*)


(** dispersion relations n(\[Lambda]) **)
(** source: Li et al. "Refractive indices of liquid crystals for display applications" J. Disp. Technol. (2005): Eq. 2 + Table VII **)
nLCeWvl[\[Lambda]_]:=Module[{ae,be,ce},   (** input \[Lambda] is in \[Mu]m **)
    {ae,be,ce} = {1.6708,0.0081,0.0024}; (** in 1, um^2, um^4 **)
    ae + be/\[Lambda]^2 + ce/\[Lambda]^4
];
nLCoWvl[\[Lambda]_]:=Module[{ao,bo,co},   (** input \[Lambda] is in \[Mu]m **)
    {ao,bo,co} = {1.5139,0.0052,0.0008}; (** in 1, um^2, um^4 **)
    ao + bo/\[Lambda]^2 + co/\[Lambda]^4
];
(** Harvey et al. "Revised Formulation for the Refractive Index of Water and Steam as a Function of Wavelength, Temperature and Density. Journal of Physical and Chemical Reference Data. 27, 761\[Dash]774 (1998) **)
(** use dimless density Overscript[\[Rho], _]=1, **)
nH2OWvl[\[Lambda]Dim_]:=Module[{\[Rho],t,\[Lambda],a0,a1,a2,a3,a4,a5,a6,a7,\[Lambda]UV,\[Lambda]IR,n},   (** input \[Lambda]Dim is in \[Mu]m **)
    {\[Rho],t,\[Lambda]}={1.0,(273.15+24)/273.15,\[Lambda]Dim/0.589};  (** dimless density, temperature, wavelength **)
    {a0,a1,a2,a3,a4,a5,a6,a7,\[Lambda]UV,\[Lambda]IR} = {0.244257733,9.74634476*^-3,-3.73234996*^-3,2.68678472*^-4,1.58920570*^-3,2.45934259*^-3,0.900704920,-1.66626219*^-2,0.2292020,5.432937}; (** dimless **)
    nVals=n/.NSolve[(n^2-1)/(n^2+2)==a0 + a1 \[Rho] + a2 t + a3 \[Lambda]^2 t + a4/\[Lambda]^2 + a5/(\[Lambda]^2-\[Lambda]UV^2) + a6/(\[Lambda]^2-\[Lambda]IR^2) + a7 \[Rho]^2 \[And] n>0,n,Reals];
	If[Length[nVals]>1,Print["Warning: More than one value found for refractive index!"];];
	nVals[[1]]
];


(* ::Subsubsection:: *)
(*plot wavelength dependence of refractive index*)


(* ::Input:: *)
(*(nLCeWvl[#*10^-3]-nLCoWvl[#*10^-3])&/@{440,550,660}*)
(*Plot[nLCeWvl[\[Lambda]*10^-3]-nLCoWvl[\[Lambda]*10^-3],{\[Lambda],400,700},Frame->True,FrameLabel->{"\[Lambda] (nm)","\[CapitalDelta]n=\!\(\*SubscriptBox[\(n\), \(e\)]\)-\!\(\*SubscriptBox[\(n\), \(o\)]\)"},Epilog->{Red,Dashed,HalfLine[{{#,0},{#,1}}]&/@{440,550,660}}]*)


(* ::Input:: *)
(*(nLCeWvl[#*10^-3])&/@{440,550,660}*)
(*(nLCoWvl[#*10^-3])&/@{440,550,660}*)
(*(nH2OWvl[#*10^-3])&/@{440,550,660}*)
(*Show[*)
(*Plot[{nLCeWvl[\[Lambda]*10^-3],nLCoWvl[\[Lambda]*10^-3]},{\[Lambda],400,700},PlotLegends->{"\!\(\*SubscriptBox[\(n\), \(LC, e\)]\)","\!\(\*SubscriptBox[\(n\), \(LC, o\)]\)"}],*)
(*Plot[nH2OWvl[\[Lambda]*10^-3],{\[Lambda],400,700},PlotLegends->{"\!\(\*SubscriptBox[\(n\), \(\*SubscriptBox[\(H\), \(2\)] O\)]\)"},PlotStyle->ColorData[97,3]],*)
(*Frame->True,Epilog->{Red,Dashed,HalfLine[{{#,0},{#,1}}]&/@{440,550,660}},*)
(*PlotRange->All,Axes->False,FrameStyle->Directive[Black]*)
(*]*)


(* ::Subsection::Closed:: *)
(*renderLCs*)


(* ::Subsubsection:: *)
(*code*)


renderLCs[solution_,{tiltAngleTheta_,tiltAnglePhi_}]:=Module[{rotaMatrix,rs,n\[Theta]s,n\[Phi]s,r,n\[Theta],n\[Phi],originalPoints,rotatedPoints,vectors,rotatedPoint,x,y,z,orientation,lcColor,liquidCrystals},
	rotaMatrix=RotationMatrix[{{0,0,1},er[tiltAngleTheta,tiltAnglePhi]}];
	
	rs={0.96,0.7,0.4,0.1};
	n\[Theta]s={7,5,3,2};
	n\[Phi]s={15,15,15,3};
	originalPoints=Flatten[MapThread[
		({r,n\[Theta],n\[Phi]}={#1,#2,#3};
		r Flatten[
			Table[
				N[{#[[1]],#[[2]],Cos[\[Theta]]}&/@(Sin[\[Theta]]CirclePoints[Round[n\[Phi] Sin[\[Theta]]+5]])]
			,{\[Theta],Subdivide[0.15,0.95 \[Pi]/2,n\[Theta]-1]}]
		,1])&
	,{rs,n\[Theta]s,n\[Phi]s}],1];
	
	rotatedPoints=rotaMatrix . #&/@originalPoints;
	vectors=Table[
		rotatedPoint=rotaMatrix . p;
		{x,y,z}=p;
		orientation=rotaMatrix . Normalize[{u[x,y,z],v[x,y,z],w[x,y,z]}/.solution];
		(rotatedPoint+# 0.05orientation)&/@{-1,1}
	,{p,originalPoints}];
	
	lcColor=RGBColor[Rational[71,85],Rational[42,85],Rational[42,85]];
	liquidCrystals=Graphics3D[{Opacity[0.8],lcColor,EdgeForm[None],Evaluate[Cylinder[#,0.02]&/@vectors]}]
]


(* ::Subsection::Closed:: *)
(*setDefaultParameters*)


setDefaultParameters[]:=Module[{},
	rDropletDim=3.5*^-6; (** m **)
	(*\[Lambda]=660.0*^-9;       (** m, red **)*)
	\[Lambda]=550.0*^-9;         (** m, green **)
	(*\[Lambda]=440.0*^-9;         (** m, blue **)*)
	noAngle=-100;
	na=0.25; (** numerical aperture **)
	magnification=1.0;
];


(* ::Subsection:: *)
(*setRefractionIndices*)


setRefractionIndices[\[Lambda]_,\[CapitalDelta]nLC0_:-1]:=Module[{},  (** input in m **)
	nOil=1.282; (** no dispersion relation data available **)
	nWater=1.333(*nH2OWvl[10^6 \[Lambda]]*); (** it's essentially constant 1.34 to 1.33 **)
	nLCe=nLCeWvl[10^6 \[Lambda]];
	nLCo=nLCoWvl[10^6 \[Lambda]];
	\[CapitalDelta]nLC=If[\[CapitalDelta]nLC0<0,nLCe-nLCo,\[CapitalDelta]nLC0];
];


(* ::Chapter::Closed:: *)
(*set physical parameters*)


setDefaultParameters[];
setRefractionIndices[\[Lambda]];


(* ::Chapter:: *)
(*compute LC alignment in upper half of droplet*)


(* ::Subsection:: *)
(*compute alignment via director (Oseen-Frank theory -> harmonic map)*)


(* ::Input:: *)
(*solution=solveHemisphereLC[];*)


(* ::Text:: *)
(*assumes spatial coordinates are rescaled by rDroplet -> droplet radius is 1*)
(*Note: Constant length constraint is not enforced sufficiently with Lagrange multiplier! -> use Q-tensor theory! (see Intro to q-tensor theory)*)


(* ::Text:: *)
(*source: W. Wang, L. Zhang, P. Zhang, Modelling and computation of liquid crystals. Acta Numerica. 30, 765\[Dash]851 (2021)*)


(* ::Subsubsection::Closed:: *)
(*test boundary conditions*)


(* ::Input:: *)
(*bcSphere=(IdentityMatrix[3]-TensorProduct[{x,y,z},{x,y,z}]) . {0,0,-1};   (** non-normalized leads to weird behavior **)*)
(*bcSphere={x z,y z,-Sqrt[1-z^2]};  (** non-normalized **)*)
(*bcSphere=Normalize[(IdentityMatrix[3]-TensorProduct[{x,y,z},{x,y,z}]) . {0,0,1}];*)
(*bcSphere={(x z)/Sqrt[x^2+y^2],(y z)/Sqrt[x^2+y^2],-Sqrt[1-z^2]};*)
(*bcSphere=If[x^2+y^2>0.1,{(x z)/Sqrt[x^2+y^2],(y z)/Sqrt[x^2+y^2],-Sqrt[1-z^2]},{1,0,0}];*)
(*bcSphere=If[x^2+y^2>0,{(x z)/Sqrt[x^2+y^2],(y z)/Sqrt[x^2+y^2],-Sqrt[1-z^2]},{(x z)/10^(-10.0),(y z)/10^(-10.0),-Sqrt[1-z^2]}];*)


(* ::Subsubsection:: *)
(*code*)


solveHemisphereLC[\[Theta]0_:0.0,\[Phi]0_:0.0]:=Module[{eps,\[Lambda]Lagrange,solution},
	Clear[u,v,w,x,y,z];
	halfBall=ImplicitRegion[x^2+y^2+z^2<=1\[And]z>=0,{x,y,z}];
	eps=10.0^-9;
	bcSphere={(x z Cos[\[Theta]0]+Sin[\[Theta]0] (-x y Cos[\[Phi]0]+(-1-eps+x^2) Sin[\[Phi]0]))/Sqrt[1+eps-(z Cos[\[Theta]0]+Sin[\[Theta]0] (-y Cos[\[Phi]0]+x Sin[\[Phi]0]))^2],(y z Cos[\[Theta]0]+Sin[\[Theta]0] ((1+eps-y^2) Cos[\[Phi]0]+x y Sin[\[Phi]0]))/Sqrt[1+eps-(z Cos[\[Theta]0]+Sin[\[Theta]0] (-y Cos[\[Phi]0]+x Sin[\[Phi]0]))^2],((-1-eps+z^2) Cos[\[Theta]0]+z Sin[\[Theta]0] (-y Cos[\[Phi]0]+x Sin[\[Phi]0]))/Sqrt[1+eps-(z Cos[\[Theta]0]+Sin[\[Theta]0] (-y Cos[\[Phi]0]+x Sin[\[Phi]0]))^2]};
	
	\[Lambda]Lagrange=0.0;  (** Problem: Fails if norm constraint is even slightly enforced! **)
	solution=NDSolve[{
	Laplacian[u[x,y,z],{x,y,z}]==-(\!\(\*SuperscriptBox[\(w\), 
TagBox[
RowBox[{"(", 
RowBox[{"0", ",", "0", ",", "1"}], ")"}],
Derivative],
MultilineFunction->None]\)[x,y,z]^2+\!\(\*SuperscriptBox[\(v\), 
TagBox[
RowBox[{"(", 
RowBox[{"0", ",", "1", ",", "0"}], ")"}],
Derivative],
MultilineFunction->None]\)[x,y,z]^2+\!\(\*SuperscriptBox[\(u\), 
TagBox[
RowBox[{"(", 
RowBox[{"1", ",", "0", ",", "0"}], ")"}],
Derivative],
MultilineFunction->None]\)[x,y,z]^2)u[x,y,z]+2 \[Lambda]Lagrange (u[x,y,z] \!\(\*SuperscriptBox[\(u\), 
TagBox[
RowBox[{"(", 
RowBox[{"1", ",", "0", ",", "0"}], ")"}],
Derivative],
MultilineFunction->None]\)[x,y,z]+v[x,y,z] \!\(\*SuperscriptBox[\(v\), 
TagBox[
RowBox[{"(", 
RowBox[{"1", ",", "0", ",", "0"}], ")"}],
Derivative],
MultilineFunction->None]\)[x,y,z]+w[x,y,z] \!\(\*SuperscriptBox[\(w\), 
TagBox[
RowBox[{"(", 
RowBox[{"1", ",", "0", ",", "0"}], ")"}],
Derivative],
MultilineFunction->None]\)[x,y,z]),
	Laplacian[v[x,y,z],{x,y,z}]==-(\!\(\*SuperscriptBox[\(w\), 
TagBox[
RowBox[{"(", 
RowBox[{"0", ",", "0", ",", "1"}], ")"}],
Derivative],
MultilineFunction->None]\)[x,y,z]^2+\!\(\*SuperscriptBox[\(v\), 
TagBox[
RowBox[{"(", 
RowBox[{"0", ",", "1", ",", "0"}], ")"}],
Derivative],
MultilineFunction->None]\)[x,y,z]^2+\!\(\*SuperscriptBox[\(u\), 
TagBox[
RowBox[{"(", 
RowBox[{"1", ",", "0", ",", "0"}], ")"}],
Derivative],
MultilineFunction->None]\)[x,y,z]^2)v[x,y,z]+2 \[Lambda]Lagrange (u[x,y,z] \!\(\*SuperscriptBox[\(u\), 
TagBox[
RowBox[{"(", 
RowBox[{"0", ",", "1", ",", "0"}], ")"}],
Derivative],
MultilineFunction->None]\)[x,y,z]+v[x,y,z] \!\(\*SuperscriptBox[\(v\), 
TagBox[
RowBox[{"(", 
RowBox[{"0", ",", "1", ",", "0"}], ")"}],
Derivative],
MultilineFunction->None]\)[x,y,z]+w[x,y,z] \!\(\*SuperscriptBox[\(w\), 
TagBox[
RowBox[{"(", 
RowBox[{"0", ",", "1", ",", "0"}], ")"}],
Derivative],
MultilineFunction->None]\)[x,y,z]),
	Laplacian[w[x,y,z],{x,y,z}]==-(\!\(\*SuperscriptBox[\(w\), 
TagBox[
RowBox[{"(", 
RowBox[{"0", ",", "0", ",", "1"}], ")"}],
Derivative],
MultilineFunction->None]\)[x,y,z]^2+\!\(\*SuperscriptBox[\(v\), 
TagBox[
RowBox[{"(", 
RowBox[{"0", ",", "1", ",", "0"}], ")"}],
Derivative],
MultilineFunction->None]\)[x,y,z]^2+\!\(\*SuperscriptBox[\(u\), 
TagBox[
RowBox[{"(", 
RowBox[{"1", ",", "0", ",", "0"}], ")"}],
Derivative],
MultilineFunction->None]\)[x,y,z]^2)w[x,y,z]+2 \[Lambda]Lagrange (u[x,y,z] \!\(\*SuperscriptBox[\(u\), 
TagBox[
RowBox[{"(", 
RowBox[{"0", ",", "0", ",", "1"}], ")"}],
Derivative],
MultilineFunction->None]\)[x,y,z]+v[x,y,z] \!\(\*SuperscriptBox[\(v\), 
TagBox[
RowBox[{"(", 
RowBox[{"0", ",", "0", ",", "1"}], ")"}],
Derivative],
MultilineFunction->None]\)[x,y,z]+w[x,y,z] \!\(\*SuperscriptBox[\(w\), 
TagBox[
RowBox[{"(", 
RowBox[{"0", ",", "0", ",", "1"}], ")"}],
Derivative],
MultilineFunction->None]\)[x,y,z]),
	DirichletCondition[Thread[{u[x,y,z],v[x,y,z],w[x,y,z]}==bcSphere],x^2+y^2+z^2==1\[And]z>0],
	DirichletCondition[Thread[{u[x,y,z],v[x,y,z],w[x,y,z]}=={0,0,-1}],z==0],
	DirichletCondition[Thread[{u[x,y,z],v[x,y,z],w[x,y,z]}=={0,0,-1}],x==0\[And]y==0]
	},{u,v,w},{x,y,z}\[Element]halfBall][[1]];
	
	solution
]


(* ::Subsubsection::Closed:: *)
(*code - quarter halfball using symmetry (fails)*)


(* ::Text:: *)
(*Neumann-0 boundary conditions are used implicitly if none are specified explicitly*)


(* ::Input:: *)
(*Clear[u,v,w,x,y,z];*)
(*halfBallQuarter=ImplicitRegion[x^2+y^2+z^2<=1\[And]x>=0\[And]y>=0\[And]z>=0,{x,y,z}];*)
(*bcSphere={(x z)/Sqrt[x^2+y^2],(y z)/Sqrt[x^2+y^2],-Sqrt[1-z^2]};*)
(*solutionSymmetry=NDSolve[{*)
(*Laplacian[u[x,y,z],{x,y,z}]==-(\!\(\*SuperscriptBox[\(w\), *)
(*TagBox[*)
(*RowBox[{"(", *)
(*RowBox[{"0", ",", "0", ",", "1"}], ")"}],*)
(*Derivative],*)
(*MultilineFunction->None]\)[x,y,z]^2+\!\(\*SuperscriptBox[\(v\), *)
(*TagBox[*)
(*RowBox[{"(", *)
(*RowBox[{"0", ",", "1", ",", "0"}], ")"}],*)
(*Derivative],*)
(*MultilineFunction->None]\)[x,y,z]^2+\!\(\*SuperscriptBox[\(u\), *)
(*TagBox[*)
(*RowBox[{"(", *)
(*RowBox[{"1", ",", "0", ",", "0"}], ")"}],*)
(*Derivative],*)
(*MultilineFunction->None]\)[x,y,z]^2)u[x,y,z],*)
(*Laplacian[v[x,y,z],{x,y,z}]==-(\!\(\*SuperscriptBox[\(w\), *)
(*TagBox[*)
(*RowBox[{"(", *)
(*RowBox[{"0", ",", "0", ",", "1"}], ")"}],*)
(*Derivative],*)
(*MultilineFunction->None]\)[x,y,z]^2+\!\(\*SuperscriptBox[\(v\), *)
(*TagBox[*)
(*RowBox[{"(", *)
(*RowBox[{"0", ",", "1", ",", "0"}], ")"}],*)
(*Derivative],*)
(*MultilineFunction->None]\)[x,y,z]^2+\!\(\*SuperscriptBox[\(u\), *)
(*TagBox[*)
(*RowBox[{"(", *)
(*RowBox[{"1", ",", "0", ",", "0"}], ")"}],*)
(*Derivative],*)
(*MultilineFunction->None]\)[x,y,z]^2)v[x,y,z],*)
(*Laplacian[w[x,y,z],{x,y,z}]==-(\!\(\*SuperscriptBox[\(w\), *)
(*TagBox[*)
(*RowBox[{"(", *)
(*RowBox[{"0", ",", "0", ",", "1"}], ")"}],*)
(*Derivative],*)
(*MultilineFunction->None]\)[x,y,z]^2+\!\(\*SuperscriptBox[\(v\), *)
(*TagBox[*)
(*RowBox[{"(", *)
(*RowBox[{"0", ",", "1", ",", "0"}], ")"}],*)
(*Derivative],*)
(*MultilineFunction->None]\)[x,y,z]^2+\!\(\*SuperscriptBox[\(u\), *)
(*TagBox[*)
(*RowBox[{"(", *)
(*RowBox[{"1", ",", "0", ",", "0"}], ")"}],*)
(*Derivative],*)
(*MultilineFunction->None]\)[x,y,z]^2)w[x,y,z],*)
(*DirichletCondition[Thread[{u[x,y,z],v[x,y,z],w[x,y,z]}==bcSphere],x^2+y^2+z^2==1\[And]x>0\[And]y>0\[And]z>0],*)
(*DirichletCondition[Thread[{u[x,y,z],v[x,y,z],w[x,y,z]}=={0,0,-1}],z==0],DirichletCondition[Thread[{u[x,y,z],v[x,y,z],w[x,y,z]}=={0,0,-1}],x==0\[And]y==0]*)
(*},{u,v,w},{x,y,z}\[Element]halfBallQuarter][[1]];*)


(* ::Subsubsection::Closed:: *)
(*derivation - constraint force*)


(* ::Input:: *)
(*Grad[\[Lambda]Lagrange(u[x,y,z]^2+v[x,y,z]^2+w[x,y,z]^2-1),{x,y,z}]//Simplify*)


(* ::Subsubsection::Closed:: *)
(*derivation*)


(* ::Input:: *)
(*Tr[Grad[{u[x,y,z],v[x,y,z],w[x,y,z]},{x,y,z}]\[Transpose]Grad[{u[x,y,z],v[x,y,z],w[x,y,z]},{x,y,z}]]*)


(* ::Text:: *)
(*derivation of vector tangential to sphere surface (Subscript[e, \[Theta]]) in cartesian coordinates*)


(* ::Input:: *)
(*Clear[x,y,z];*)
(*e\[Theta]Cart0=TransformedField["Spherical"->"Cartesian",{0,1,0}, {r,\[Theta],\[CurlyPhi]}->{x,y,z}]*)
(*e\[Theta]Cart=((Numerator[e\[Theta]Cart0]/.#)/(Denominator[e\[Theta]Cart0]/.#))&@{Sqrt[x^2+y^2+z^2]->1,Sqrt[x^2+y^2]->Sqrt[1-z^2]}*)
(*norm=Sqrt[# . #]&@e\[Theta]Cart/.{z^2->1-x^2-y^2,z^4->(1-x^2-y^2)^2}//Simplify*)
(*norm=(ComplexExpand[Norm[#]]&@(e\[Theta]Cart/.{z^2->1-x^2-y^2,z^4->(1-x^2-y^2)^2})//Simplify)/.{x^2+y^2+z^2->1}*)


(* ::Input:: *)
(*Clear[r,x,y,z];*)
(*e\[Theta]Cart0=TransformedField["Spherical"->"Cartesian",{0,1,0}, {r,\[Theta],\[CurlyPhi]}->{x,y,z}];*)
(*((Numerator[e\[Theta]Cart0]/.#)/(Denominator[e\[Theta]Cart0]/.#))&@{Sqrt[x^2+y^2+z^2]->1,Sqrt[x^2+y^2]->Sqrt[1-z^2]}*)
(*e\[Theta]Cart[{x_,y_,z_}]:={(x z)/Sqrt[1-z^2],(y z)/Sqrt[1-z^2],-Sqrt[1-z^2]}*)
(*rotaMatrix=RotationMatrix[\[Pi]/2+\[Phi],{0,0,1}] . RotationMatrix[\[Theta],{1,0,0}];*)
(*rotaMatrixInv=RotationMatrix[-\[Theta],{1,0,0}] . RotationMatrix[-(\[Pi]/2+\[Phi]),{0,0,1}];*)
(*rotaMatrix . e\[Theta]Cart[rotaMatrixInv . {x,y,z}]//Simplify*)


(* ::Input:: *)
(*Clear[x,y,z];*)
(*rotaMatrix . {x,y,z};*)
(*e\[Theta]Cart[{x,y,z}];*)
(*{\[Theta]0,\[Phi]0,eps}={1/3 \[Pi]/2,\[Pi]/2,10.0^-9};*)
(*Simplify[e\[Theta]Cart[rotaMatrix . {x,y,z}]]/.{\[Theta]->\[Theta]0,\[Phi]->\[Phi]0,1->1+eps};*)
(*Simplify[rotaMatrix . e\[Theta]Cart[{x,y,z}]]/.{\[Theta]->\[Theta]0,\[Phi]->\[Phi]0,1->1+eps};*)
(*bc3d=Simplify[rotaMatrix . e\[Theta]Cart[rotaMatrixInv . {x,y,z}]/.{\[Theta]->\[Theta]0,\[Phi]->\[Phi]0,1->1+eps}];*)
(*bc3d[[{1,2}]]/.z->Sqrt[1-x^2-y^2]*)
(*VectorPlot[bc3d[[{1,2}]]/.z->Sqrt[1-x^2-y^2],{x,-1,1},{y,-1,1},RegionFunction->Function[{x,y,vx,vy,n},x^2+y^2<=0.9]]*)
(*points=Table[{Cos[\[Phi]],Sin[\[Phi]]},{\[Phi],Subdivide[0,\[Pi],10-1]}];*)
(*directions=Map[( *)
(*{xx,zz}=#;*)
(*Normalize[bc3d[[{1,3}]]/.{x->xx,y->0,z->zz}]*)
(*)&,points];*)
(*vectorPoints=MapThread[{#1-0.1 0.5#2,#1+0.1 0.5#2}&,{points,directions}];*)
(*Graphics[{Black,Arrow/@vectorPoints}]*)
(**)
(*points=Table[{Cos[\[Phi]],Sin[\[Phi]]},{\[Phi],Subdivide[0,\[Pi],10-1]}];*)
(*directions=Map[( *)
(*{yy,zz}=#;*)
(*Normalize[bc3d[[{2,3}]]/.{x->0,y->yy,z->zz}]*)
(*)&,points];*)
(*vectorPoints=MapThread[{#1-0.1 0.5#2,#1+0.1 0.5#2}&,{points,directions}];*)
(*Graphics[{Black,Arrow/@vectorPoints}]*)


(* ::Input:: *)
(*Manipulate[*)
(*Graphics3D[{Sphere[],Arrow[Line[{er[\[Theta],\[Phi]],er[\[Theta],\[Phi]]+e\[Theta][\[Theta],\[Phi]]}]]}]*)
(*,{\[Theta],0,\[Pi],\[Pi]/10.}*)
(*,{\[Phi],0,2\[Pi],2\[Pi]/10.}]*)


(* ::Subsection::Closed:: *)
(*compute alignment via tilt angle*)


(* ::Text:: *)
(*get tilt angle as a linear change from top to bottom boundary*)
(*-> curvature = 0 -> \[CapitalDelta]\[Theta]=0*)


(* ::Subsubsection::Closed:: *)
(*test boundary conditions*)


(* ::Input:: *)
(*bcSphere=(IdentityMatrix[3]-TensorProduct[{x,y,z},{x,y,z}]) . {0,0,-1};   (** non-normalized leads to weird behavior **)*)
(*bcSphere={x z,y z,-Sqrt[1-z^2]};  (** non-normalized **)*)
(*bcSphere=Normalize[(IdentityMatrix[3]-TensorProduct[{x,y,z},{x,y,z}]) . {0,0,1}];*)
(*bcSphere={(x z)/Sqrt[x^2+y^2],(y z)/Sqrt[x^2+y^2],-Sqrt[1-z^2]};*)
(*bcSphere=If[x^2+y^2>0.1,{(x z)/Sqrt[x^2+y^2],(y z)/Sqrt[x^2+y^2],-Sqrt[1-z^2]},{1,0,0}];*)
(*bcSphere=If[x^2+y^2>0,{(x z)/Sqrt[x^2+y^2],(y z)/Sqrt[x^2+y^2],-Sqrt[1-z^2]},{(x z)/10^(-10.0),(y z)/10^(-10.0),-Sqrt[1-z^2]}];*)


(* ::Subsubsection::Closed:: *)
(*code 1d*)


(* ::Code:: *)
(*Clear[r];*)
(*height[r_]:=Sqrt[1-r^2];*)
(*bcTop=ArcTan[D[height[r],r]];*)
(*radius=0.5;*)
(*\[Theta]Anchor=0;*)
(*bcBottom=\[Theta]Anchor;*)
(*solution\[Theta]=NDSolve[{*)
(*D[\[Theta][z],{z,2}]==0,*)
(*\[Theta][height[radius]]==bcTop/.r->radius,*)
(*\[Theta][0]==bcBottom*)
(*},\[Theta],{z,0,height[radius]}][[1]]*)


(* ::Subsubsection::Closed:: *)
(*code 3d*)


(* ::Code:: *)
(*Clear[x,y,z];*)
(*halfBall=ImplicitRegion[x^2+y^2+z^2<=1\[And]z>=0,{x,y,z}];*)
(*height[r_]:=Sqrt[1-r^2];*)
(*bcTop[r_]:=ArcTan[D[height[r],r]];*)
(*radius=0.5;*)
(*\[Theta]Anchor=0;*)
(*bcBottom=\[Theta]Anchor;*)
(*solution\[Theta]=NDSolve[{*)
(*Laplacian[\[Theta][x,y,z],{x,y,z}]==0,*)
(*\[Theta][x,y,z]==bcTop/.r->radius,*)
(*\[Theta][x,y,z]==bcBottom*)
(*},\[Theta],{x,y,z}\[Element]halfBall][[1]]*)


(* ::Subsubsection::Closed:: *)
(*code - quarter halfball using symmetry (fails)*)


(* ::Text:: *)
(*Neumann-0 boundary conditions are used implicitly if none are specified explicitly*)


(* ::Input:: *)
(*Clear[u,v,w,x,y,z];*)
(*halfBallQuarter=ImplicitRegion[x^2+y^2+z^2<=1\[And]x>=0\[And]y>=0\[And]z>=0,{x,y,z}];*)
(*bcSphere={(x z)/Sqrt[x^2+y^2],(y z)/Sqrt[x^2+y^2],-Sqrt[1-z^2]};*)
(*solutionSymmetry=NDSolve[{*)
(*Laplacian[u[x,y,z],{x,y,z}]==-(\!\(\*SuperscriptBox[\(w\), *)
(*TagBox[*)
(*RowBox[{"(", *)
(*RowBox[{"0", ",", "0", ",", "1"}], ")"}],*)
(*Derivative],*)
(*MultilineFunction->None]\)[x,y,z]^2+\!\(\*SuperscriptBox[\(v\), *)
(*TagBox[*)
(*RowBox[{"(", *)
(*RowBox[{"0", ",", "1", ",", "0"}], ")"}],*)
(*Derivative],*)
(*MultilineFunction->None]\)[x,y,z]^2+\!\(\*SuperscriptBox[\(u\), *)
(*TagBox[*)
(*RowBox[{"(", *)
(*RowBox[{"1", ",", "0", ",", "0"}], ")"}],*)
(*Derivative],*)
(*MultilineFunction->None]\)[x,y,z]^2)u[x,y,z],*)
(*Laplacian[v[x,y,z],{x,y,z}]==-(\!\(\*SuperscriptBox[\(w\), *)
(*TagBox[*)
(*RowBox[{"(", *)
(*RowBox[{"0", ",", "0", ",", "1"}], ")"}],*)
(*Derivative],*)
(*MultilineFunction->None]\)[x,y,z]^2+\!\(\*SuperscriptBox[\(v\), *)
(*TagBox[*)
(*RowBox[{"(", *)
(*RowBox[{"0", ",", "1", ",", "0"}], ")"}],*)
(*Derivative],*)
(*MultilineFunction->None]\)[x,y,z]^2+\!\(\*SuperscriptBox[\(u\), *)
(*TagBox[*)
(*RowBox[{"(", *)
(*RowBox[{"1", ",", "0", ",", "0"}], ")"}],*)
(*Derivative],*)
(*MultilineFunction->None]\)[x,y,z]^2)v[x,y,z],*)
(*Laplacian[w[x,y,z],{x,y,z}]==-(\!\(\*SuperscriptBox[\(w\), *)
(*TagBox[*)
(*RowBox[{"(", *)
(*RowBox[{"0", ",", "0", ",", "1"}], ")"}],*)
(*Derivative],*)
(*MultilineFunction->None]\)[x,y,z]^2+\!\(\*SuperscriptBox[\(v\), *)
(*TagBox[*)
(*RowBox[{"(", *)
(*RowBox[{"0", ",", "1", ",", "0"}], ")"}],*)
(*Derivative],*)
(*MultilineFunction->None]\)[x,y,z]^2+\!\(\*SuperscriptBox[\(u\), *)
(*TagBox[*)
(*RowBox[{"(", *)
(*RowBox[{"1", ",", "0", ",", "0"}], ")"}],*)
(*Derivative],*)
(*MultilineFunction->None]\)[x,y,z]^2)w[x,y,z],*)
(*DirichletCondition[Thread[{u[x,y,z],v[x,y,z],w[x,y,z]}==bcSphere],x^2+y^2+z^2==1\[And]x>0\[And]y>0\[And]z>0],*)
(*DirichletCondition[Thread[{u[x,y,z],v[x,y,z],w[x,y,z]}=={0,0,-1}],z==0],DirichletCondition[Thread[{u[x,y,z],v[x,y,z],w[x,y,z]}=={0,0,-1}],x==0\[And]y==0]*)
(*},{u,v,w},{x,y,z}\[Element]halfBallQuarter][[1]];*)


(* ::Subsubsection::Closed:: *)
(*derivation*)


(* ::Input:: *)
(*Tr[Grad[{u[x,y,z],v[x,y,z],w[x,y,z]},{x,y,z}]\[Transpose]Grad[{u[x,y,z],v[x,y,z],w[x,y,z]},{x,y,z}]]*)


(* ::Text:: *)
(*derivation of vector tangential to sphere surface (Subscript[e, \[Theta]]) in cartesian coordinates*)


(* ::Input:: *)
(*Clear[x,y,z];*)
(*e\[Theta]Cart0=TransformedField["Spherical"->"Cartesian",{0,1,0}, {r,\[Theta],\[CurlyPhi]}->{x,y,z}]*)
(*e\[Theta]Cart=((Numerator[e\[Theta]Cart0]/.#)/(Denominator[e\[Theta]Cart0]/.#))&@{Sqrt[x^2+y^2+z^2]->1,Sqrt[x^2+y^2]->Sqrt[1-z^2]}*)
(*norm=Sqrt[# . #]&@e\[Theta]Cart/.{z^2->1-x^2-y^2,z^4->(1-x^2-y^2)^2}//Simplify*)
(*norm=(ComplexExpand[Norm[#]]&@(e\[Theta]Cart/.{z^2->1-x^2-y^2,z^4->(1-x^2-y^2)^2})//Simplify)/.{x^2+y^2+z^2->1}*)


(* ::Input:: *)
(*Manipulate[*)
(*Graphics3D[{Sphere[],Arrow[Line[{er[\[Theta],\[Phi]],er[\[Theta],\[Phi]]+e\[Theta][\[Theta],\[Phi]]}]]}]*)
(*,{\[Theta],0,\[Pi],\[Pi]/10.}*)
(*,{\[Phi],0,2\[Pi],2\[Pi]/10.}]*)


(* ::Subsection::Closed:: *)
(*set alignment via tilt angle, analytical*)


rd=rDropletDim/rDropletDim;  (** rescaled droplet radius **)

\[Theta]rLinear[r_]:=\[Theta]Center + r (\[Theta]Radius-\[Theta]Center)/rd
\[Alpha]=10/rd;
\[Theta]rExp[r_]:=\[Theta]Center-((\[Theta]Center-\[Theta]Radius) Exp[-\[Alpha] (rd-r)]+\[Theta]Radius)
Plot[{\[Theta]rLinear[r],\[Theta]rExp[r]},{r,0,rd},Frame->True,FrameLabel->{"radius \!\(\*
StyleBox[\"r\",\nFontSlant->\"Italic\"]\)","tilt \!\(\*
StyleBox[\"\[Theta]\",\nFontSlant->\"Italic\"]\)"},FrameStyle->Directive[Black,18]];

h[r_]:=Sqrt[rd^2-r^2];
\[Beta][r_]:=10.0/h[r];
\[Theta]zLinear[z_,r_]:=Module[{\[Theta]Bottom,\[Theta]Top},
	\[Theta]Bottom=\[Theta]rExp[r];
	\[Theta]Top=\[Pi]/2-ArcTan[-h'[r]];
	(\[Theta]Top-\[Theta]Bottom)z/h[r]+\[Theta]Bottom
]
\[Theta]zExp[z_,r_]:=Module[{\[Theta]Bottom,\[Theta]Top},
	\[Theta]Bottom=\[Theta]rExp[r];
	\[Theta]Top=\[Pi]/2-ArcTan[-h'[r]];
	(\[Theta]Top-\[Theta]Bottom)Exp[-\[Beta][r](h[r]-z)]+\[Theta]Bottom
]

With[{r=0.1rd},Plot[{\[Theta]zLinear[z,r],\[Theta]zExp[z,r]},{z,0,h[r]},Frame->True,FrameLabel->{"position \!\(\*
StyleBox[\"z\",\nFontSlant->\"Italic\"]\)","tilt \!\(\*
StyleBox[\"\[Theta]\",\nFontSlant->\"Italic\"]\)"},FrameStyle->Directive[Black,18]]];

directorAnalytical[{x_,y_,z_}]:=Module[{r},
	r=Norm[{x,y}];
	{\[Theta]Bottom,\[Theta]Top}={\[Theta]rExp[r],\[Pi]/2-ArcTan[-h'[r]]};
	er[\[Theta]zExp[z,r],If[x==0\[And]y==0,0.0,\[Pi]+ArcTan[x,y]]]
]


(* ::Subsubsection::Closed:: *)
(*2d*)


(* ::Input:: *)
(*xzVector;*)


(* ::Input:: *)
(*xz\[Theta]=Table[*)
(*r=Abs[x];*)
(*{\[Theta]Bottom,\[Theta]Top}={\[Theta]rExp[r],\[Pi]/2-ArcTan[-h'[r]]};*)
(*Table[{x,z,\[Theta]zExp[z,Abs@x]If[x<0,1,-1]},{z,Subdivide[0,h[r],10]}]*)
(*,{x,Subdivide[-0.99rd,0.99rd,20]}];*)
(*xzVector={#[[{1,2}]],{Sin[#[[3]]],Cos[#[[3]]]}}&/@Flatten[xz\[Theta],1];*)
(**)
(*ListVectorPlot[xzVector,AspectRatio->1/2,Frame->True,Axes->False,Epilog->{Gray,Opacity[0.5],Circle[{0,0},rd]},VectorPoints->Flatten[xz\[Theta][[All,All,{1,2}]],1],VectorMarkers->"Segment",FrameLabel->{"x (\[Mu]m)","z (\[Mu]m)"}]*)


(* ::Subsubsection::Closed:: *)
(*3d*)


(* ::Input:: *)
(*Flatten[xyzDirector[[All,All,All]],2]*)


(* ::Input:: *)
(*xyzDirector=Table[*)
(*r=Norm[{x,y}];*)
(*If[r<rd,*)
(*{\[Theta]Bottom,\[Theta]Top}={\[Theta]rExp[r],\[Pi]/2-ArcTan[-h'[r]]};*)
(*Print[{\[Theta]Bottom,\[Theta]Top}];*)
(*(*Table[{{x,y,z},er[\[Theta]zExp[z],If[x==0\[And]y==0,0.0,ArcTan[x,y]]]},{z,Subdivide[0,h[r],10]}]*)*)
(*Table[\[Theta]zExp[z],{z,Subdivide[0,h[r],10]}]*)
(*,Sequence@@{}]*)
(*,{x,Subdivide[-0.99rd,0.99rd,10][[{10}]]}*)
(*,{y,{0}(*Subdivide[-0.99rd,0.99rd,10]*)}]*)


(* ::Input:: *)
(*xyzDirector=Table[*)
(*r=Norm[{x,y}];*)
(*If[r<rd,*)
(*{\[Theta]Bottom,\[Theta]Top}={\[Theta]rExp[r],\[Pi]/2-ArcTan[-h'[r]]};*)
(*Table[{{x,y,z},er[\[Theta]zExp[z,r],If[x==0\[And]y==0,0.0,\[Pi]+ArcTan[x,y]]]},{z,Subdivide[0,h[r],10-1]}]*)
(*,Sequence@@{}]*)
(*(*,{x,(*{0}*)Subdivide[-0.99rd,0.99rd,10-1]}*)
(*,{y,{0}(*Subdivide[-0.99rd,0.99rd,10-1]*)}];*)*)
(*(*,{x,{0}}*)
(*,{y,Subdivide[-0.99rd,0.99rd,10-1]}];*)*)
(*,{x,Subdivide[-0.99rd,0.99rd,10-1]}*)
(*,{y,Subdivide[-0.99rd,0.99rd,10-1]}];*)
(*lcVectors=Map[{#[[1]]- 0.05#[[2]],#[[1]]+ 0.05#[[2]]}&,Flatten[xyzDirector[[All,All,All]],2]];*)


(* ::Input:: *)
(*xyzDirector=Table[*)
(*{x,y}={Cos[\[Phi]],Sin[\[Phi]]};*)
(*Print[{\[Phi],N@ArcTan[x,y],-N@ArcTan[x,y]}];*)
(*r=Norm[{x,y}];*)
(*If[r<rd,*)
(*{\[Theta]Bottom,\[Theta]Top}={\[Theta]rExp[r],\[Pi]/2-ArcTan[-h'[r]]};*)
(*Table[{{x,y,z},er[\[Theta]zExp[z],If[x==0\[And]y==0,0.0,\[Pi]+ArcTan[x,y]]]},{z,Subdivide[0,h[r],10-1]}]*)
(*,Sequence@@{}]*)
(*(*,{x,(*{0}*)Subdivide[-0.99rd,0.99rd,10-1]}*)
(*,{y,{0}(*Subdivide[-0.99rd,0.99rd,10-1]*)}];*)*)
(*,{\[Phi],Subdivide[0,2\[Pi],10][[;;-2]]}];*)
(*lcVectors=Map[{#[[1]]- 0.05#[[2]],#[[1]]+ 0.05#[[2]]}&,Flatten[xyzDirector[[All,All]],1]];*)


(* ::Input:: *)
(*plot=Graphics3D[{Evaluate[CapsuleShape[#,0.03]&/@lcVectors],Opacity[0.2],Sphere[{0,0,0},rd]},PlotRange->{All,All,{0,All}},Boxed->False,ViewVertical->{-0.11293301187983383`,0.0,1.0},ViewPoint->{-2.4149686489954596`,-1.5612429571636024`,1.7833807369925514`}]*)


(* ::Subsection::Closed:: *)
(*save director*)


(* ::Subsubsection::Closed:: *)
(*extended boundary function*)


(* ::Input:: *)
(*bcSphereTop[{x_,y_,z_}]:={(x z)/Sqrt[x^2+y^2+0.000000001],(y z)/Sqrt[x^2+y^2+0.000000001],-(Norm[{x,y,z}]^2-z^2)/Sqrt[x^2+y^2+0.000000001]}*)


(* ::Input:: *)
(*r=1.5;*)
(*y=0;*)
(*VectorPlot[bcSphereTop[{x,y,z}][[{1,3}]],{x,-2,2},{z,0,2},AspectRatio->1/2,Prolog->{Gray,Circle[{0,0},r]}]*)


(* ::Input:: *)
(*bcBottom={0,0,-1}*)
(*bcHemisphere[{x_,y_,z_}]:={(x z)/Sqrt[x^2+y^2+0.000000001],(y z)/Sqrt[x^2+y^2+0.000000001],-(1^2-z^2)/Sqrt[x^2+y^2+0.000000001]}*)


(* ::Input:: *)
(*directorApproximation[{x_,y_,z_}]:=Normalize[{x,y,z}]*)


(* ::Input:: *)
(*directorApproximation[{x_,y_,z_}]:=If[Norm[{x,y}]<1,Normalize[bcBottom+z/Sqrt[1-Norm[{x,y}]^2+0.000000001](bcHemisphere[{x,y,z}]-bcBottom)],{0,0,0}]*)


(* ::Input:: *)
(*r=1.0;*)
(*y=0;*)
(*VectorPlot[directorApproximation[{x,y,z}][[{1,3}]],{x,-1.1,1.1},{z,0,1.1},AspectRatio->1/2,Prolog->{Gray,Circle[{0,0},r]}]*)


(* ::Input:: *)
(*directorApproximation[{0.9,y,0.1}]//Norm*)


(* ::Input:: *)
(*r=1.0;*)
(*y=0;*)
(*StreamPlot[directorApproximation[{x,y,z}][[{1,3}]],{x,-1.1,1.1},{z,0,1.1},AspectRatio->1/2,Prolog->{Gray,Circle[{0,0},r]}]*)


(* ::Subsubsection::Closed:: *)
(*.h5 for meep*)


(* ::Input:: *)
(*(** regular bulk beyond boundaries **)*)
(*radius=1.0;*)
(*xVals=Subdivide[-1.05radius,1.05radius,70];*)
(*yVals=Subdivide[-1.05radius,1.05radius,70];*)
(*zVals=Subdivide[-0.05radius,1.05radius,36];*)
(*samplePoints=Flatten[Table[{x,y,z},{x,xVals},{y,yVals},{z,zVals}],2];*)
(*Length@samplePoints*)
(*(*Graphics3D[{Red,Point[samplePoints],(*Black,Point[boundarySamplePoints],*)Opacity[0.2],Sphere[]}]*)*)


(* ::Input:: *)
(*directorVals//Dimensions*)


(* ::Input:: *)
(*directorVals=*)
(*Table[If[z>0,If[Norm[{x,y,z}]<0.99radius,director[{x,y,z}],bcSphereTop[{x,y,z}]],{0,0,-1}],{x,xVals},{y,yVals},{z,zVals}];(*If[#\[LeftDoubleBracket]3\[RightDoubleBracket]>0,If[Norm[#]<0.99radius,director[#],bcSphereTop[#]],{0,0,-1}]&/@samplePoints;*)*)
(*Export[FileNameJoin[{NotebookDirectory[],"positionDirectorRegularGrid.hdf5"}],{xVals,yVals,zVals,directorVals},{"Datasets",{"xVals","yVals","zVals","directorVals"}}]*)


(* ::Subsubsection::Closed:: *)
(*.h5 for meep (old)*)


(* ::Input:: *)
(*(** regular bulk and boundary **)*)
(*radius=1.0;*)
(*bulkSamplePoints=Flatten[Table[If[Norm[{x,y,z}]<=1.05radius\[And]z>-0.05,{x,y,z},Nothing],{x,Subdivide[-1.05radius,1.05radius,70]},*)
(*{y,Subdivide[-1.05radius,1.05radius,70]},*)
(*{z,Subdivide[-1.05radius,1.05radius,70]}]*)
(*,2];*)
(*boundarySamplePoints=If[#[[3]]>0.01,#,Nothing]&/@(radius SpherePoints[2000]);*)
(*samplePoints=bulkSamplePoints;(*Join[bulkSamplePoints,boundarySamplePoints];*)*)
(*Length/@{bulkSamplePoints,boundarySamplePoints}*)
(*Graphics3D[{Red,Point[bulkSamplePoints],(*Black,Point[boundarySamplePoints],*)Opacity[0.2],Sphere[]}]*)


(* ::Input:: *)
(*(** random bulk and regular boundary **)*)
(*radius=0.995;*)
(*bulkSamplePoints=RandomPoint[ImplicitRegion[x^2+y^2+z^2<radius^2\[And]z>=0,{x,y,z}],3000];*)
(*boundarySamplePoints=If[#[[3]]>0.01,#,Nothing]&/@(radius SpherePoints[2000]);*)
(*samplePoints=Join[bulkSamplePoints,boundarySamplePoints];*)
(*Length/@{bulkSamplePoints,boundarySamplePoints}*)
(*Graphics3D[{Red,Point[bulkSamplePoints],Black,Point[boundarySamplePoints]}]*)


(* ::Input:: *)
(*(** random bulk and boundary: fails **)*)
(*radius=0.995;*)
(*bulkSamplePoints=RandomPoint[ImplicitRegion[x^2+y^2+z^2<radius^2\[And]z>=0,{x,y,z}],3000];*)
(*boundarySamplePoints=RandomPoint[ImplicitRegion[x^2+y^2+z^2==radius^2\[And]z>=0,{{x,-radius,radius},{y,-radius,radius},{z,0,radius}}],3000]*)
(*samplePoints=Join[bulkSamplePoints,boundarySamplePoints];*)
(*Graphics3D[{Blue,Point[samplePoints]}]*)


(* ::Subsubsection::Closed:: *)
(*.dat*)


(* ::Input:: *)
(*lookupTablePositionDirector=Flatten[Table[{{x,y,z},If[z>0,If[Norm[{x,y,z}]<1,director[x,y,z],{0,0,0}],{0,0,0}]},{x,0,1,0.2},{y,-1,1,0.2},{z,0.01,1,0.2}],2];*)
(*Export[FileNameJoin[{NotebookDirectory[],"directorLookupTable.dat"}],lookupTablePositionDirector];*)


(* ::Subsection::Closed:: *)
(*visualize director*)


(* ::Text:: *)
(*also see figures.wl Chapter Figure S1: Optics*)


(* ::Subsubsection::Closed:: *)
(*symmetry: cut through volume - continuous field*)


(* ::Input:: *)
(*yCut=0.0;*)
(*quarterDisk=ImplicitRegion[x^2+z^2<=1\[And]x>=0\[And]z>=0,{x,z}];*)
(*Clear[x,z];*)
(*DensityPlot[Evaluate[u[x,yCut,z]/.solutionSymmetry],{x,z}\[Element]quarterDisk,PlotLegends->Placed[Automatic,Right]]*)
(*DensityPlot[Evaluate[v[x,yCut,z]/.solutionSymmetry],{x,z}\[Element]quarterDisk,PlotLegends->Placed[Automatic,Right]]*)
(*DensityPlot[Evaluate[w[x,yCut,z]/.solutionSymmetry],{x,z}\[Element]quarterDisk,PlotLegends->Placed[Automatic,Right]]*)
(*DensityPlot[Evaluate[ArcTan[u[x,yCut,z],w[x,yCut,z]]/.solutionSymmetry],{x,z}\[Element]quarterDisk,AspectRatio->1,PlotLegends->Placed[Automatic,Right]]*)
(*DensityPlot[Evaluate[ArcTan[v[yCut,x,z],w[yCut,x,z]]/.solutionSymmetry],{x,z}\[Element]quarterDisk,AspectRatio->1,PlotLegends->Placed[Automatic,Right]]*)
(*StreamPlot[Evaluate[{u[x,yCut,z],w[x,yCut,z]}/.solutionSymmetry],{x,z}\[Element]quarterDisk,AspectRatio->1,PlotLegends->Placed[Automatic,Right],PlotRange->All,StreamPoints->Fine]*)


(* ::Subsubsection::Closed:: *)
(*cut through volume - continuous field*)


(* ::Input:: *)
(*yCut=0.0;*)
(*halfDisk=ImplicitRegion[x^2+z^2<=1\[And]z>=0,{x,z}];*)
(*StreamPlot[Evaluate[{u[x,yCut,z],w[x,yCut,z]}/.solution],{x,z}\[Element]halfDisk,AspectRatio->1/2,PlotLegends->Placed[Automatic,Right],PlotRange->All,StreamPoints->Fine]*)


(* ::Input:: *)
(*yCut=0.0;*)
(*halfDisk=ImplicitRegion[x^2+z^2<=1\[And]z>=0,{x,z}];*)
(*Clear[x,z];*)
(*DensityPlot[Evaluate[u[x,yCut,z]/.solution],{x,z}\[Element]halfDisk,AspectRatio->1/2,PlotLegends->Placed[Automatic,Right]]*)
(*DensityPlot[Evaluate[v[x,yCut,z]/.solution],{x,z}\[Element]halfDisk,AspectRatio->1/2,PlotLegends->Placed[Automatic,Right]]*)
(*DensityPlot[Evaluate[w[x,yCut,z]/.solution],{x,z}\[Element]halfDisk,AspectRatio->1/2,PlotLegends->Placed[Automatic,Right]]*)
(*DensityPlot[Evaluate[ArcTan[u[x,yCut,z],w[x,yCut,z]]/.solution],{x,z}\[Element]halfDisk,AspectRatio->1/2,PlotLegends->Placed[Automatic,Right]]*)
(*DensityPlot[Evaluate[ArcTan[v[yCut,x,z],w[yCut,x,z]]/.solution],{x,z}\[Element]halfDisk,AspectRatio->1/2,PlotLegends->Placed[Automatic,Right]]*)
(*StreamPlot[Evaluate[{u[x,yCut,z],w[x,yCut,z]}/.solution],{x,z}\[Element]halfDisk,AspectRatio->1/2,PlotLegends->Placed[Automatic,Right],PlotRange->All,StreamPoints->Fine]*)


(* ::Subsubsection::Closed:: *)
(*cut through volume - LC*)


(* ::Input:: *)
(*{nxSamples,nzSamples}={20,5};*)


(* ::Input:: *)
(*circleSamplePoints=If[Sqrt[#[[1]]^2+#[[2]]^2]<0.99^2\[And]#[[2]]>=0,#,Sequence@@{}]&/@Flatten[Outer[{##}&,N@Subdivide[-1,1,nxSamples-1],N@Subdivide[-1,1,nzSamples-1]],1];*)
(*circleSamplePoints//Length*)
(*ListPlot[circleSamplePoints,PlotRange->{{-1,1},{-1,1}},AspectRatio->1,Frame->True]*)


(* ::Input:: *)
(*circleSamplePoints=Flatten[Table[Table[{x,z}*)
(*,{z,Subdivide[0,Sqrt[0.99^2-x^2],nzSamples-1]}]*)
(*,{x,N@Subdivide[-1,1,nxSamples-1]}],1];*)
(*circleSamplePoints//Length*)
(*Show[*)
(*ListPlot[circleSamplePoints,PlotRange->{{-1,1},{-1,1}},AspectRatio->1,Frame->True,Epilog->{Black,Circle[]}],*)
(*Plot[Sqrt[0.99^2-x^2],{x,-1,1}]*)
(*]*)


(* ::Subsubsection::Closed:: *)
(*cut through rotated volume - LC*)


(* ::Input:: *)
(*{nxSamples,nzSamples}={20,10};*)


(* ::Input:: *)
(*originalPoints=Flatten[Table[*)
(*zMax=Sqrt[effectiveRadius^2-x^2];*)
(*zMin=1-effectiveRadius;*)
(*If[zMax>zMin,Table[{x,yCut,z},{z,N@Subdivide[zMin,zMax,nzSamples-1]}],Sequence@@{}]*)
(*,{x,N@Subdivide[-effectiveRadius,effectiveRadius,nxSamples-1]}],1];*)
(*rotatedPoints=rotaMatrix . #&/@originalPoints;*)
(*Graphics3D[{Blue,Point[originalPoints],Red,Point[rotatedPoints]}]*)


(* ::Input:: *)
(*(*tiltAnglesTheta=N@{\[Pi]/4};*)
(*tiltAnglesPhi=N@{\[Pi]/4};*)*)
(*tiltAnglesTheta=N@{0};*)
(*tiltAnglesPhi=N@{0};*)
(*tiltAngleTheta=tiltAnglesTheta[[1]];*)
(*tiltAnglePhi=tiltAnglesPhi[[1]];*)
(*rotaMatrix=RotationMatrix[{{0,0,1},tiltDirector[tiltAngleTheta,tiltAnglePhi]}];*)
(*rotaMatrixInv=RotationMatrix[{tiltDirector[tiltAngleTheta,tiltAnglePhi],{0,0,1}}];*)
(**)
(*effectiveRadius=0.99;*)
(*yCut=0.0;*)
(*AbsoluteTiming[*)
(*{nxSamples,nzSamples}={20,10};*)
(*originalPoints=Flatten[Table[*)
(*zMax=Sqrt[effectiveRadius^2-x^2];*)
(*zMin=1-effectiveRadius;*)
(*If[zMax>zMin,Table[{x,yCut,z},{z,N@Subdivide[zMin,zMax,nzSamples-1]}],Sequence@@{}]*)
(*,{x,N@Subdivide[-effectiveRadius,effectiveRadius,nxSamples-1]}],1];*)
(*rotatedPoints=rotaMatrix . #&/@originalPoints;*)
(*vectors=Table[*)
(*rotatedPoint=rotaMatrix . p;*)
(*{x,y,z}=p;*)
(*orientation=rotaMatrix . Normalize[{u[x,y,z],v[x,y,z],w[x,y,z]}/.solution];*)
(*(rotatedPoint+# 0.05orientation)&/@{-1,1}*)
(*(*(p+# 0.05director[p])&/@{-1,1}*)*)
(*,{p,originalPoints}];*)
(*][[1]]*)
(**)
(*plot=Graphics3D[*)
(*{Darker@Red,Evaluate[CapsuleShape[#,0.02]&/@vectors],Opacity[0.2],Lighter@Red,Sphere[]},PlotRange->{All,All,{-1,1}},Boxed->False,ViewVertical->{0.045639271885001866`,-0.37241498404788465`,0.9269434376047012`},ViewPoint->{0.4026446942044853`,-3.355677274521836`,0.16524914367911303`}]*)


(* ::Subsubsection::Closed:: *)
(*volume (3d LC alignment)*)


(* ::Input:: *)
(*AbsoluteTiming[*)
(*{nrMax,nz,n\[Phi]}={10,4,15};*)
(*(*{nrMax,nz,n\[Phi]}={10,16,75};*)*)
(*points=Flatten[*)
(*Table[Table[*)
(*angle=ArcSin[z];*)
(*N[{#[[1]],#[[2]],Sin[angle]}&/@(0.95r Cos[angle]CirclePoints[Round[n\[Phi] r Cos[angle]+5]])],{r,Subdivide[0.0,0.99,Round[nrMax(1-z)]]}],{z,Subdivide[0.1,0.95,nz-1]}],2];*)
(*vectors=Table[*)
(*{x,y,z}=p;*)
(*orientation=Normalize[{u[x,y,z],v[x,y,z],w[x,y,z]}/.solution];*)
(*(p+# 0.05orientation)&/@{-1,1}*)
(*(*(p+# 0.05director[p])&/@{-1,1}*)*)
(*,{p,points}];*)
(*(*vectors=Map[( *)
(*orientation=director[#];*)
(*{#- 0.05orientation,#+ 0.05orientation}*)
(*)&,points];*)*)
(*][[1]]*)


(* ::Input:: *)
(*plot=Graphics3D[{(*Darker@Red,*)EdgeForm[None],Evaluate[Cylinder[#,0.03]&/@vectors],Opacity[0.2],Sphere[]},PlotRange->{All,All,{0,All}},Boxed->False,ViewVertical->{-0.11293301187983383`,0.0,1.0},ViewPoint->{-2.4149686489954596`,-1.5612429571636024`,1.7833807369925514`}]*)


(* ::Input:: *)
(*Export[FileNameJoin[{NotebookDirectory[],"directorDistribution.png"}],plot]*)


(* ::Chapter:: *)
(*Optics / Jones calculus*)


(* ::Subsection:: *)
(*Jones calculus*)


(* ::Subsubsection::Closed:: *)
(*Jones vectors*)


polarizationStateHorizontal={1,0};
polarizationStateVertical={0,1};
polarizationStateDiagonal=1.0/Sqrt[2.0]{1,1};
polarizationVector={1,0,0};
orthogonalPolarizationVector={0,1,0};
waveVector={0,0,-1};


(* ::Subsubsection:: *)
(*Jones matrices*)


setJonesCalculusDiscretization[rDropletDim_]:=Module[{nxJones,nyJones,nzJones,dzJones,xVals,yVals,zVals,d},
	{nxJones,nyJones,nzJones}=101{1,1,1};    (* 101{1,1,2} *)
	xVals=Subdivide[-1.0,1.0,nxJones-1];
	yVals=Subdivide[-1.0,1.0,nyJones-1];
	zVals=Subdivide[-1.0,1.0,nzJones-1];
	dzJones=(zVals[[2]]-zVals[[1]]);
	d=rDropletDim dzJones;                   (** scales distance element from 2 to 2 rDropletDim **)
	{xVals,yVals,zVals,nxJones,nyJones,nzJones,d}
];


{xVals,yVals,zVals,nxJones,nyJones,nzJones,d}=setJonesCalculusDiscretization[rDropletDim];

\[Delta]LC[\[CapitalDelta]n_]:= 2\[Pi]/\[Lambda] d \[CapitalDelta]n;     (** \[CapitalDelta]n = ne - no; liquid crystal birefringent uniaxial medium **)
\[Delta]LC[\[CapitalDelta]n_,d_]:= 2\[Pi]/\[Lambda] d \[CapitalDelta]n;  (**  **)
\[Delta]Iso[n_,d_]:= 2\[Pi]/\[Lambda] d n;   (** optically isotropic medium **)
(*jonesMatrixLC[\[Theta]_,\[Delta]e_]:=Exp[\[ImaginaryI] \[Delta]Iso[nLCo]]{{Cos[\[Theta]]^2+Sin[\[Theta]]^2 #,Sin[\[Theta]]Cos[\[Theta]](#-1.0)},{Sin[\[Theta]]Cos[\[Theta]](#-1.0),Sin[\[Theta]]^2+Cos[\[Theta]]^2 #}}&@Exp[\[ImaginaryI] \[Delta]e];*)
(** usage: jonesMatrixLC[\[Theta],\[Delta]LC[\[CapitalDelta]nLC],\[Delta]Iso[nLCo]] **)
jonesMatrixLC[\[Theta]_,\[Delta]e_,\[Delta]o_]:=Exp[I \[Delta]o]{{Cos[\[Theta]]^2+Sin[\[Theta]]^2 #,Sin[\[Theta]]Cos[\[Theta]](#-1.0)},{Sin[\[Theta]]Cos[\[Theta]](#-1.0),Sin[\[Theta]]^2+Cos[\[Theta]]^2 #}}&@Exp[I \[Delta]e];

jonesMatrixFC[d_]:=Exp[I \[Delta]Iso[nOil,d]]IdentityMatrix[2];
jonesMatrixH2O[d_]:=Exp[I \[Delta]Iso[nWater,d]]IdentityMatrix[2];

(** usage: jonesMatrixIsotropic[\[Delta]Uni[\[CapitalDelta]n,d]] **)
jonesMatrixIsotropic[\[Delta]_]:=Exp[I \[Delta]]IdentityMatrix[2];
jonesMatrixHorizontalPolarizer={{1,0},{0,0}};
jonesMatrixVerticalPolarizer={{0,0},{0,1}};
jonesMatrixPolarizer[\[Theta]_]:=N@RotationMatrix[\[Theta]] . jonesMatrixHorizontalPolarizer . RotationMatrix[-\[Theta]]


jonesCalculusPropagator[xVals_,yVals_,directorAnglesRefIndices3d_,nLCo_,d_,rDropletDim_,nxJones_,nyJones_,nzJones_,\[Lambda]_,compiledQ_:True]:=Module[{x,y,totalJonesMatrix,directorAnglesRefIndices,singleJonesMatrices,directorAngle,refIndex,beforeDropletState,afterDropletState},
If[compiledQ,
	jonesCalculusPropagatorC[xVals,yVals,nxJones,nyJones,nzJones,directorAnglesRefIndices3d,\[Lambda],rDropletDim,d,nLCo]
,
Table[
{x,y}={xVals[[xIndex]],yVals[[yIndex]]};   (** dimless **) 
totalJonesMatrix=If[Norm[{x,y}]<1,
(** ray intersection with droplet **)
(** TODO: precompute ray intersections to only use Jones matrix in LC **)
directorAnglesRefIndices=Reverse@directorAnglesRefIndices3d[[xIndex,yIndex]];
singleJonesMatrices=Map[(
{directorAngle,refIndex}=#;
If[-\[Pi]<=directorAngle<=\[Pi],
jonesMatrixLC[directorAngle,\[Delta]LC[refIndex-nLCo],\[Delta]Iso[nLCo,d]],
jonesMatrixIsotropic[\[Delta]Iso[refIndex,d]]
]
)&,directorAnglesRefIndices];
Dot@@singleJonesMatrices   (** output **)
,
(** no ray intersection with droplet: only water **)
jonesMatrixH2O[2.0rDropletDim]   (** output **)
];
beforeDropletState=jonesMatrixHorizontalPolarizer . {1.0,1.0} initialPhaseFactor[x rDropletDim,y rDropletDim];
afterDropletState=totalJonesMatrix . beforeDropletState;
jonesMatrixVerticalPolarizer . afterDropletState
(*jonesMatrixPolarizer[0.95\[Pi]/2].afterDropletState*)
,{xIndex,nxJones}
,{yIndex,nyJones}
]
]
];


(* ::Subsubsection::Closed:: *)
(*compiled code*)


jonesMatrixIsotropicC=Compile[{{d,_Real,0},{\[Lambda],_Real,0},{n,_Real,0}},
	Exp[I 2.0\[Pi]/\[Lambda] d n ]{{1.0,0.0},{0.0,1.0}}
	,CompilationTarget->"C",RuntimeOptions->"Speed",Parallelization->True
];

jonesMatrixLCC=Compile[{{\[Theta],_Real,0},{\[Lambda],_Real,0},{d,_Real,0},{nLCo,_Real,0},{nLCeMod,_Real,0}},
Module[{eio,eie},
	eie=Exp[I 2.0\[Pi]/\[Lambda] d (nLCeMod-nLCo)];
	eio=Exp[I 2.0\[Pi]/\[Lambda] d nLCo];
	eio{{Cos[\[Theta]]^2+Sin[\[Theta]]^2 eie,Sin[\[Theta]]Cos[\[Theta]](eie-1.0)},{Sin[\[Theta]]Cos[\[Theta]](eie-1.0),Sin[\[Theta]]^2+Cos[\[Theta]]^2 eie}}
]
,CompilationTarget->"C",RuntimeOptions->"Speed",Parallelization->True
];

With[{jonesMatrixHorizontalPolarizer=jonesMatrixHorizontalPolarizer,jonesMatrixVerticalPolarizer=jonesMatrixVerticalPolarizer,nWater=nWater},
jonesCalculusPropagatorC=Compile[{{xVals,_Real,1},{yVals,_Real,1},{nxJones,_Real,0},{nyJones,_Real,0},{nzJones,_Real,0},{directorAnglesRefIndices3d,_Real,4},{\[Lambda],_Real,0},{rDropletDim,_Real,0},{d,_Real,0},{nLCo,_Real,0}},
Module[{x,y,beforeDropletState,afterDropletState,totalJonesMatrix(*=Table[0.0+0.0\[ImaginaryI],{2,2}]*),
directorAnglesRefIndices,directorAngle,refIndex,singleJonesMatrices},
Table[
	{x,y}={xVals[[xIndex]],yVals[[yIndex]]};   (** dimless **) 
	If[x*x+y*y<1,
		(** ray intersection with droplet **)
	
		directorAnglesRefIndices=Reverse[directorAnglesRefIndices3d[[xIndex,yIndex]]];
		totalJonesMatrix={{1.0+0.0I,0.0},{0.0,1.0}};
		Do[
			{directorAngle,refIndex}=directorAnglesRefIndices[[i]];
			totalJonesMatrix=If[-\[Pi]<=directorAngle<=\[Pi],
				jonesMatrixLCC[directorAngle,\[Lambda],d,nLCo,refIndex],
				jonesMatrixIsotropicC[d,\[Lambda],refIndex]
			] . totalJonesMatrix;
		,{i,nzJones}];
	,
		(** no ray intersection with droplet: only water **)
		totalJonesMatrix=jonesMatrixIsotropicC[2.0rDropletDim,\[Lambda],nWater]   (** output **)
	];
	beforeDropletState=jonesMatrixHorizontalPolarizer . {1.0,1.0};
	afterDropletState=totalJonesMatrix . beforeDropletState;
	jonesMatrixVerticalPolarizer . afterDropletState
	,{xIndex,nxJones}
	,{yIndex,nyJones}
]
]
,CompilationTarget->"C",RuntimeOptions->"Speed",Parallelization->True
,CompilationOptions->{"InlineCompiledFunctions"->True,"InlineExternalDefinitions"->True}
];
];


(* ::Subsubsection::Closed:: *)
(*initial phase distribution (planar/spherical wave)*)


constantPhase[x_,y_]:=0;                                         (** planar wave in z-direction **)
sphericalWavePhase[x_,y_,z0_]:=2\[Pi]/\[Lambda] Sqrt[x^2+y^2+(0-z0)^2];     (** \[Phi]; spherical wave at z=0; z=e^(\[ImaginaryI]|k||r|); using m **)
initialPhaseFactor[x_,y_]:=Exp[I #]&@constantPhase[x,y](*sphericalWavePhase[x,y,rDropletDim]*)(*sphericalWavePhase[x,y,-1rDropletDim]*); (**  **)


(* ::Input:: *)
(*(** visualization **)*)
(*phase=Table[Arg@initialPhaseFactor[x rDropletDim,y rDropletDim],{x,Subdivide[-1,1,100]},{y,Subdivide[-1,1,100]}];*)
(*ListDensityPlot[phase,PlotLabel->"wrapped phase"]*)


(* ::Subsubsection::Closed:: *)
(*derivations*)


(* ::Input:: *)
(*Clear[\[Delta]]*)
(*Inverse[RotationMatrix[\[Theta]]] . {{1,0},{0,Exp[I \[Delta]]}} . RotationMatrix[\[Theta]]//Simplify*)
(*RotationMatrix[\[Alpha]] . IdentityMatrix[2] . RotationMatrix[-\[Alpha]]//Simplify*)
(*RotationMatrix[\[Alpha]] . DiagonalMatrix[{a,b}] . RotationMatrix[-\[Alpha]]//Simplify*)


(* ::Input:: *)
(*{{Exp[I Subscript[\[Delta], o]],0},{0,Exp[I Subscript[\[Delta], e]]}}//MatrixForm*)


(* ::Input:: *)
(*{{1,0},{0,Exp[I \[Delta]]}}//MatrixForm*)


(* ::Input:: *)
(*{{Cos[\[Theta]]^2+E^(I \[Delta]) Sin[\[Theta]]^2,(-1+E^(I \[Delta])) Cos[\[Theta]] Sin[\[Theta]]},{(-1+E^(I \[Delta])) Cos[\[Theta]] Sin[\[Theta]],E^(I \[Delta]) Cos[\[Theta]]^2+Sin[\[Theta]]^2}}//MatrixForm*)
(*{{Cos[\[Theta]]^2+Sin[\[Theta]]^2 #,Sin[\[Theta]]Cos[\[Theta]](#-1.0)},{Sin[\[Theta]]Cos[\[Theta]](#-1.0),Sin[\[Theta]]^2+Cos[\[Theta]]^2 #}}&@Exp[I \[Delta]]//MatrixForm*)


(* ::Subsubsection::Closed:: *)
(*comparison Jones matrix vs. optical path length**)


(* ::Input:: *)
(*jonesMatrixLC[0\[Degree],\[Pi]/2,0.0]*)
(*Arg@Diagonal[jonesMatrixLC[0\[Degree],\[Pi]/2,0]]*)
(*jonesMatrixLC[45\[Degree],\[Pi]/2,0.0]*)


(* ::Input:: *)
(*polarizationStateHorizontal;*)
(*Arg[polarizationStateHorizontal]*)
(*eFieldQWP=jonesMatrixLC[45\[Degree],\[Pi]/2,0.0] . polarizationStateHorizontal*)
(*Arg[eFieldQWP]*)


(* ::Input:: *)
(*{\[Phi]x0,\[Phi]y0}={0.0,0.0};*)
(*\[CapitalDelta]nLC=nLCe-nLCo;*)
(*phaseRetardation=2\[Pi]/\[Lambda] \[CapitalDelta]nLC d;*)
(*{nLCeQWP,nLCoQWP}={\[CapitalDelta]nLC (\[Pi]/2)/phaseRetardation,0.0}+nLCo*)
(*{nLCx,nLCy}=calcLCrefIndicesXY[lcDirector,xPolarizationVector,yPolarizationVector,{nLCeQWP,nLCoQWP}];*)
(*{\[Phi]h,\[Phi]v}={\[Phi]x0,\[Phi]y0}+2\[Pi]/\[Lambda] d{nLCx,nLCy}*)


(* ::Input:: *)
(*Graphics[{{EdgeForm[Black],White,Ellipsoid[{0,0},{nLCe,nLCo}]},{Black,Arrow[{{0,0},{nLCe,0}}],Arrow[{{0,0},{0,nLCo}}]}},*)
(*Frame->True,*)
(*FrameStyle->Directive[Black,AbsoluteThickness[1.5],18],*)
(*FrameLabel->{"extraordinary index \!\(\*SubscriptBox[\(n\), \(e\)]\)","ordinary index \!\(\*SubscriptBox[\(n\), \(o\)]\)"},*)
(*AspectRatio->nLCo/nLCe*)
(*]*)


(* ::Subsubsection::Closed:: *)
(*quarter wave plate test*)


(* ::Text:: *)
(*source: Osterman, PhD thesis, (2005), ch. 5, p. 29*)
(*phase retardation \[CapitalGamma]=\[Pi]/2*)
(*action: Subscript[H, QWP](45\[Degree])|Horizontal\[RightAngleBracket] = |Circular,+\[RightAngleBracket] *)
(*|Horizontal\[RightAngleBracket] = (1,0)*)
(*|Circular,+\[RightAngleBracket] = {1/2+I/2,-(1/2)+I/2}*)
(*Subscript[H, QWP](0\[Degree]) = {{Exp[-I \[Pi]/4],0},{0,Exp[I \[Pi]/4]}}*)


(* ::Input:: *)
(*(** first component is delayed, 2^nd is sped up **)*)
(*jonesMatrixQWP[\[Phi]_]:=RotationMatrix[-\[Phi]] . {{Exp[-I \[Pi]/4],0},{0,Exp[I \[Pi]/4]}} . RotationMatrix[\[Phi]];*)
(*jonesMatrixQWPnoOffset[\[Phi]_]:=RotationMatrix[-\[Phi]] . {{1,0},{0,Exp[I \[Pi]/2]}} . RotationMatrix[\[Phi]];*)
(*jtest=jonesMatrixQWPnoOffset[45\[Degree]] . polarizationStateHorizontal//Simplify*)
(*Mod[Differences[Arg/@jtest],2\[Pi]][[1]]*)
(**)
(*{nLCeTest,nLCoTest}={nLCe,nLCo};*)
(*\[CapitalDelta]nLC=nLCe-nLCo*)
(*phaseRetardation=2\[Pi]/\[Lambda] \[CapitalDelta]nLC d;*)
(*\[CapitalDelta]nLCqwp=\[CapitalDelta]nLC (\[Pi]/2)/phaseRetardation(** set optical anisotropy \[CapitalDelta]n such that LC element becomes a QWP **)*)
(*\[Pi]/2==2\[Pi]/\[Lambda] \[CapitalDelta]nLCqwp d;*)
(*jonesVecCircTest=jonesMatrixLC[45\[Degree],\[Delta]LC[\[CapitalDelta]nLCqwp],0(*\[Delta]Iso[nLCoTest]*)] . polarizationStateHorizontal*)
(*Mod[Differences[Arg/@jonesVecCircTest],2\[Pi]][[1]]==\[Pi]/2*)


(* ::Input:: *)
(*Solve[\[Pi]/2==2\[Pi]/\[Lambda] \[CapitalDelta]nLC dd,dd][[1]]*)


(* ::Input:: *)
(*argForm[n_]:=(Row@{Abs[n],Superscript["\[ExponentialE]","\[ImaginaryI] \[Pi] "<>ToString[Arg[n]/\[Pi]]]})*)


(* ::Input:: *)
(*Map[argForm[#]&,jonesMatrixLC[0\[Degree],\[Delta]LC[\[CapitalDelta]nLCqwp],(*\[Delta]Iso=*)-\[Pi]/4],{2}]*)


(* ::Input:: *)
(*Map[{Abs[#],Arg[#]/\[Pi]}&,jonesMatrixLC[0\[Degree],\[Delta]LC[\[CapitalDelta]nLCqwp],(*\[Delta]Iso=*)-\[Pi]/4],{2}]*)


(* ::Subsection::Closed:: *)
(*calculate refractive index n / optical length 2d*)


(* ::Input:: *)
(*Clear[x,y,z,director,nLC];*)
(*rayDirector={0,0,-1};*)
(*polarizationVector={1,0,0};*)
(*director[x_,y_,z_]=Normalize[{u[x,y,z],v[x,y,z],w[x,y,z]}/.solution];*)
(*(*nLCBad[x_,y_,z_]=Norm[surfaceProjection[-{0,0,1},director[x,y,z]]];*)*)
(*(*nLCBad2[cos\[Theta]_]:=1.0/Sqrt[cos\[Theta]^2/nLCo^2+(1-cos\[Theta]^2)/nLCe^2];*)*)
(*nLC[director_,polarizationVector_]:=Module[{rota,rotatedP},*)
(*rota=RotationMatrix[{director,{1.0,0,0}}];  (** align director with x-axis for simpler calculation **)*)
(*rotatedP=rota . polarizationVector;*)
(*((rotatedP[[1]]/nLCe)^2+(rotatedP[[2]]/nLCo)^2+(rotatedP[[3]]/nLCo)^2)^(-0.5)*)
(*]*)


(* ::Input:: *)
(*nLC[{0,0,1},{0,1,0}]*)


(* ::Input:: *)
(*xx=0.01;*)
(*{dx,dy,dz}=0.02{1,1,1};*)
(*tiltAngle=(*0 \[Pi]*)\[Pi]/4;*)
(**)
(*refractiveIndexData2d=Table[*)
(*refIndex=If[Norm[{xx,y,z}]<=0.99,If[z Cos[tiltAngle]+y Sin[tiltAngle]>0,*)
(*rotatedDirector=RotationMatrix[-tiltAngle,{1,0,0}] . director[xx,y Cos[tiltAngle]-z Sin[tiltAngle],z Cos[tiltAngle]+y Sin[tiltAngle]];*)
(*(*directorTilt=ArcCos[{0,0,-1}.rotatedDirector];*)*)
(*(*cosDirectorTilt=rayDirector.rotatedDirector;*)*)
(*(*nLC[cosDirectorTilt]*)*)
(*nLC[rotatedDirector,polarizationVector],nOil],nWater];*)
(*{y,z,refIndex}*)
(*,{y,-1,1,dy},{z,-1,1,dz}];*)
(*Export[FileNameJoin[{NotebookDirectory[],"refractiveIndex2d_tiltAngle_"<>ToString[DecimalForm[N@tiltAngle,{3,2}]]<>".dat"}],refractiveIndexData2d];*)
(**)
(*plot=ListDensityPlot[Flatten[refractiveIndexData2d,1],InterpolationOrder->0,FrameLabel->{"x","z"},PlotRange->All,PlotLegends->BarLegend[Automatic,LegendLabel->Placed["n",Right]](*,ColorFunctionScaling->False,ColorFunction->(ColorData["M10DefaultDensityGradient"][Rescale[#,{nOil,1.7}]]&)*),PlotRangePadding->None]*)
(*Histogram[Flatten[refractiveIndexData2d[[All,All,3]]],{nOil,1.02nLCe,0.01},Frame->True,Prolog->{Gray,InfiniteLine[{{#,0},{#,1}}]&/@{nLCo,nLCe}}]*)
(**)
(*dzz=10^6 rDropletDim dz; (** integral element **)*)
(*pathLengths=Total/@refractiveIndexData2d[[All,All,3]] dzz;*)
(*plot=ListLinePlot[{Rescale[Range[Length[pathLengths]],{1,Length[pathLengths]},{-1.0,1}],pathLengths}\[Transpose],Frame->True,FrameLabel->{"x (\!\(\*SubscriptBox[\(r\), \(droplet\)]\))","optical pathlength (\[Mu]m)"},FrameStyle->Directive[Black,14,AbsoluteThickness[2]],Axes->False]*)


(* ::Subsection::Closed:: *)
(*calculate refractive indices as a function of tilt angle*)


(* ::Input:: *)
(*Clear[x,y,z,director,nLC];*)
(*director[{x_?NumericQ,y_?NumericQ,z_?NumericQ}]:=Normalize[{u[x,y,z],v[x,y,z],w[x,y,z]}/.solution];*)
(*(*refIndexLC[x_,y_,z_]=Norm[surfaceProjection[-{0,0,1},director[x,y,z]]];*)*)
(*nLC[cos\[Theta]_]:=1.0/Sqrt[cos\[Theta]^2/nLCo^2+(1-cos\[Theta]^2)/nLCe^2];*)


(* ::Input:: *)
(*xx=0.01;*)
(*{dx,dy,dz}=0.02{1,1,1};*)
(*tiltAngles=N@Subdivide[0,\[Pi]/2,20];*)
(*Print[AbsoluteTiming[*)
(*Do[*)
(*refractiveIndexData2d=Table[*)
(*refIndex=If[Norm[{xx,y,z}]<=0.99,If[z Cos[tiltAngle]+y Sin[tiltAngle]>0,*)
(*rotatedDirector=RotationMatrix[-tiltAngle,{1,0,0}] . director[xx,y Cos[tiltAngle]-z Sin[tiltAngle],z Cos[tiltAngle]+y Sin[tiltAngle]];*)
(*cosDirectorTilt={0,0,-1} . rotatedDirector;*)
(*(*directorTilt=ArcCos[{0,0,-1}.rotatedDirector];*)*)
(*nLC[cosDirectorTilt],nOil],nWater];*)
(*{y,z,refIndex}*)
(*,{y,-1,1,dy},{z,-1,1,dz}];*)
(*Export[FileNameJoin[{NotebookDirectory[],"refractiveIndex2d_tiltAngle_"<>ToString[DecimalForm[N@tiltAngle,{3,2}]]<>".mat"}],refractiveIndexData2d];*)
(*,{tiltAngle,tiltAngles}];*)
(*][[1]]];*)


(* ::Subsection::Closed:: *)
(*calculate polarization state*)


(* ::Text:: *)
(*assume that polarization angle point in x-direction is = 0 = 2\[Pi]*)


(* ::Input:: *)
(*xx=0.01;*)
(*{dx,dy,dz}=0.02{1,1,1};*)
(*{nx,ny,nz}=10{1,1,1};*)
(*xVals=Subdivide[-1.0,1.0,nx-1];*)
(*yVals=Subdivide[-1.0,1.0,ny-1];*)
(*zVals=Subdivide[-1.0,1.0,nz-1];*)
(*dz=(zVals[[2]]-zVals[[1]]);*)
(*tiltAngle=0 \[Pi](*\[Pi]/4*);*)
(**)
(*polarizationAngles=Table[*)
(*polarizationAngle=If[Norm[{xx,y,z}]<=0.99,*)
(*If[z Cos[tiltAngle]+y Sin[tiltAngle]>0,*)
(*rotatedDirector=RotationMatrix[-tiltAngle,{1,0,0}] . director[xx,y Cos[tiltAngle]-z Sin[tiltAngle],z Cos[tiltAngle]+y Sin[tiltAngle]];*)
(*polarizationAngleLC=ArcTan[rotatedDirector[[1]],rotatedDirector[[2]]];*)
(*polarizationAngleLC,-10],-10];*)
(*{y,z,polarizationAngle}*)
(*,{y,-1,1,dy},{z,-1,1,dz}];*)


(* ::Input:: *)
(*noAngle=-10;*)
(*Print[AbsoluteTiming[*)
(*polarizationAngles3d=Table[*)
(*polarizationAngle=If[Norm[{x,y,z}]<=0.99,*)
(*If[z Cos[tiltAngle]+y Sin[tiltAngle]>0,*)
(*rotatedDirector=RotationMatrix[-tiltAngle,{1,0,0}] . director[x,y Cos[tiltAngle]-z Sin[tiltAngle],z Cos[tiltAngle]+y Sin[tiltAngle]];*)
(*polarizationAngleLC=ArcTan[rotatedDirector[[1]],rotatedDirector[[2]]];*)
(*polarizationAngleLC,noAngle],noAngle];*)
(*{x,y,z,polarizationAngle}*)
(*,{x,-1,1,dx},{y,-1,1,dy},{z,-1,1,dz}];*)
(*][[1]]];*)


(* ::Input:: *)
(*xVals=N@Subdivide[-1,1,100];*)
(*yVals=N@Subdivide[-1,1,100];*)
(*zVals=N@Subdivide[-1,1,100];*)
(*noAngle=-10;*)
(*Print[AbsoluteTiming[*)
(*polarizationAngles3dParallel=ParallelTable[*)
(*polarizationAngle=If[Norm[{x,y,z}]<=0.99,*)
(*If[z Cos[tiltAngle]+y Sin[tiltAngle]>0,*)
(*rotatedDirector=RotationMatrix[-tiltAngle,{1,0,0}] . director[x,y Cos[tiltAngle]-z Sin[tiltAngle],z Cos[tiltAngle]+y Sin[tiltAngle]];*)
(*polarizationAngleLC=ArcTan[rotatedDirector[[1]],rotatedDirector[[2]]];*)
(*polarizationAngleLC,noAngle],noAngle];*)
(*{x,y,z,polarizationAngle}*)
(*,{x,xVals},{y,yVals},{z,zVals}];*)
(*][[1]]];*)


(* ::Input:: *)
(*Histogram[Flatten@Mod[#,2\[Pi]]&@polarizationAngles]*)


(* ::Input:: *)
(*plot=ListDensityPlot[Flatten[polarizationAngles,1],InterpolationOrder->0,FrameLabel->{"x","z"},PlotRange->All,PlotLegends->BarLegend[Automatic,LegendLabel->Placed["n",Right]](*,ColorFunctionScaling->False,ColorFunction->(ColorData["M10DefaultDensityGradient"][Rescale[#,{nOil,1.7}]]&)*),PlotRangePadding->None]*)


(* ::Input:: *)
(*dropletPolarizationState=Table[*)
(*startState={1,0};*)
(*filteredAngles=Select[polarizationAngles3d[[x,y,All,4]],#!=noAngle&];*)
(*singleJonesMatrices=jonesMatrix[#]&/@filteredAngles;*)
(*dropletJonesMatrix=Dot@@singleJonesMatrices;*)
(*afterDropletState=dropletJonesMatrix . startState;*)
(*afterPolarizerState=jonesMatrix[\[Pi]/2] . afterDropletState;*)
(*intensity=Abs[afterPolarizerState[[2]]]^2;*)
(*intensity*)
(*,{x,101},{y,101}];*)


(* ::Input:: *)
(*ListDensityPlot[dropletPolarizationState,InterpolationOrder->0,PlotRangePadding->None,PlotLegends->BarLegend[Automatic,LegendLabel->"intensity"],FrameLabel->{"x","y"}(*,ColorFunction->Hue*)]*)


(* ::Input:: *)
(*avgAngle=Table[*)
(*filteredAngles=Select[polarizationAngles3d[[x,y,All,4]],#!=noAngle&];*)
(*Mean[filteredAngles]*)
(*,{x,101},{y,101}];*)


(* ::Input:: *)
(*ListDensityPlot[avgAngle,InterpolationOrder->0,PlotRangePadding->None,PlotLegends->BarLegend[Automatic,LegendLabel->"avg. angle"],FrameLabel->{"x","y"},ColorFunction->Hue]*)


(* ::Subsubsection::Closed:: *)
(*test*)


(* ::Input:: *)
(*surfaceProjection[v_,nUnit_]:=(IdentityMatrix[Length[v]]-TensorProduct[#,#]&@nUnit) . v*)


(* ::Input:: *)
(*Solve[(nee Cos[\[Theta]])^2/no^2+(nee Sin[\[Theta]])^2/ne^2==1,nee]*)


(* ::Input:: *)
(*directorYData2d=Table[*)
(*refIndex=If[Norm[{xx,y,z}]<=1,If[z>0,(*Norm[director[xx,y,z]]*)(*director[xx,y,z]\[LeftDoubleBracket]3\[RightDoubleBracket]*)(*ArcCos[{0,0,1}.director[xx,y,z]]*)(*{0,0,1}.director[xx,y,z]*){0,1,0} . director[xx,y,z],nOil],nWater];*)
(*{y,z,refIndex}*)
(*,{y,-1,1,dy},{z,-1,1,dz}];*)
(*plot=ListDensityPlot[Flatten[directorYData2d,1],InterpolationOrder->0,FrameLabel->{"x","z"},PlotRange->All,PlotLegends->BarLegend[Automatic](*,ColorFunctionScaling->False,ColorFunction->(ColorData["M10DefaultDensityGradient"][Rescale[#,{1.2,2}]]&)*)];*)


(* ::Input:: *)
(*xx=0.7;*)
(*{dy,dz}={0.1,0.1};*)
(*refractiveIndexData2d=Table[*)
(*refIndex=If[Norm[{xx,y,z}]<=1,If[z>0,refIndexLC[xx,y,z],refIndexOil],refIndexWater];*)
(*{y,z,refIndex}*)
(*,{y,-1,1,dy},{z,-1,1,dz}];*)
(*Export[FileNameJoin[{NotebookDirectory[],"refractiveIndex2d.dat"}],refractiveIndexData2d];*)
(*plot=ListDensityPlot[Flatten[refractiveIndexData2d,1],InterpolationOrder->0,FrameLabel->{"x","z"},PlotRange->All,PlotLegends->Automatic]*)
(**)
(*opticalLengthData=Flatten[Table[*)
(*Flatten[{{x,y,z},If[Norm[{x,y,z}]<=1,If[z>0,refIndexLC[x,y,z],refIndexOil],refIndexWater]}]*)
(*,{x,0.001,1,0.2},{y,-1,1,0.2},{z,-1,1,0.2}],2];*)


(* ::Subsection::Closed:: *)
(*maltese cross as a function of tilt angle \[Theta]**)


(* ::Text:: *)
(*deprecated; TODO: calculate refractive indices*)


(* ::Subsubsection:: *)
(*main*)


(* ::Input:: *)
(*tiltAngles={0}(*N[Subdivide[0,\[Pi]/2,20]]*);*)
(*tiltAnglePhi=0.0(*\[Pi]/4.0*);*)
(*Print[AbsoluteTiming[*)
(*intensityTiltData=Table[*)
(*Print[AbsoluteTiming[*)
(*rotaMatrix=RotationMatrix[{{0,0,1},tiltDirector[tiltAngle,tiltAnglePhi]}];*)
(*rotaMatrixInv=RotationMatrix[{tiltDirector[tiltAngle,tiltAnglePhi],{0,0,1}}];*)
(*polarizationAngles3d=Table[*)
(*position={x,y,z};*)
(*polarizationAngle=If[Norm[position]<=0.99,*)
(*If[position . tiltDirector[tiltAngle,tiltAnglePhi]>0.0,*)
(*rotatedPosition=rotaMatrixInv . position;*)
(*rotatedDirector=rotaMatrix . (director@@rotatedPosition);*)
(*polarizationAngleLC=ArcTan[rotatedDirector[[1]],rotatedDirector[[2]]];*)
(*polarizationAngleLC,lastAngle],noAngle];*)
(*{x,y,z,polarizationAngle}*)
(*,{x,xVals},{y,yVals},{z,zVals}];*)
(*][[1]]];*)
(**)
(*dropletPolarizationState=Table[*)
(*startState={1,0};*)
(*filteredAngles=Select[polarizationAngles3d[[x,y,All,4]],#!=noAngle&];*)
(*(*singleJonesMatrices=If[#!=lastAngle,jonesMatrixLC[#-nLCo],jonesMatrixFC]&/@filteredAngles;*)*)
(*singleJonesMatrices=jonesMatrixLC[\[Delta]LC[#-nLCo],\[Delta]Iso[nLCo]]&/@filteredAngles;*)
(*dropletJonesMatrix=Dot@@singleJonesMatrices;*)
(*afterDropletState=dropletJonesMatrix . startState;*)
(*afterPolarizerState=jonesMatrixLC[\[Pi]/2] . afterDropletState;*)
(*intensity=(*Abs[afterPolarizerState\[LeftDoubleBracket]1\[RightDoubleBracket]]^2+*)Abs[afterPolarizerState[[2]]]^2;*)
(*{afterPolarizerState,intensity}*)
(*,{x,nx},{y,ny}];*)
(**)
(*{tiltAngle,Mean[Flatten@dropletPolarizationState[[2]]],dropletPolarizationState}*)
(*,{tiltAngle,tiltAngles}];*)
(*][[1]]];*)


(* ::Input:: *)
(*plot=Show[*)
(*Plot[0.35Cos[\[Theta]]^2,{\[Theta],0,\[Pi]/2},PlotStyle->LightGray,PlotLegends->{"\!\(\*SuperscriptBox[\(cos\), \(2\)]\)\[Theta]"}],*)
(*ListLinePlot[intensityTiltData[[All,{1,2}]],Mesh->All,PlotLegends->{"Jones matrix"}]*)
(*,Frame->True,FrameLabel->{"tilt angle \[Theta]","intensity"},PlotLabel->Style["\[Phi] = "<>ToString[Round[tiltAnglePhi 180/\[Pi]]]<>"\[Degree]",Black,20],FrameStyle->Directive[Black,18,AbsoluteThickness[2]],AspectRatio->1,PlotRange->{0,0.4}*)
(*]*)
(*Export[FileNameJoin[{NotebookDirectory[],"intensity_over_tilt_phi_"<>ToString[IntegerString[tiltAnglePhi 180/\[Pi]]]<>".png"}],plot]*)


(* ::Input:: *)
(*i=-1;*)
(*ListDensityPlot[intensityTiltData[[i,3,All,All,2]],InterpolationOrder->0,PlotRangePadding->None,PlotRange->All,PlotLegends->BarLegend[Automatic,LegendLabel->"intensity"],FrameLabel->{"x","y"},PlotLabel->"\[Theta] = "<>ToString[DecimalForm[intensityTiltData[[i,1]]180/\[Pi],{4,1}]]<>" \[Degree]",ColorFunction->GrayLevel,PerformanceGoal->"Speed"]*)


(* ::Input:: *)
(*counter=0;*)
(*Do[*)
(*label="\[Theta] = "<>ToString[DecimalForm[intensityTiltData[[i,1]]180/\[Pi],{4,1}]]<>"\[Degree], \[Phi] = "<>ToString[DecimalForm[tiltAnglePhi 180/\[Pi],{4,1}]]<>"\[Degree]";*)
(*plotHighRes=ListDensityPlot[intensityTiltData[[i,3]],InterpolationOrder->0,PlotRangePadding->None,PlotRange->All,PlotLegends->BarLegend[Automatic,LegendLabel->"intensity"],FrameLabel->{"x","y"},ColorFunction->GrayLevel,ImageSize->{Automatic,400}];*)
(*plotLowRes=ListDensityPlot[intensityTiltData[[i,3]],InterpolationOrder->0,PlotRangePadding->None,PlotRange->All,PlotLegends->BarLegend[Automatic,LegendLabel->"intensity"],FrameLabel->{"x","y"},ColorFunction->GrayLevel,ImageSize->{Automatic,400},PerformanceGoal->"Speed"];*)
(*plot=Grid[{{label,SpanFromLeft},{plotHighRes,plotLowRes}}];*)
(*Export[FileNameJoin[{NotebookDirectory[],"mma_polarizationMovie_"<>ToString[tiltAnglePhi 180/\[Pi]]<>"_deg","p_"<>IntegerString[++counter,10,IntegerLength[Length@intensityTiltData]+1]<>".png"}],plot];*)
(*,{i,Length[intensityTiltData]}];*)


(* ::Input:: *)
(*Manipulate[*)
(*ListDensityPlot[intensityTiltData[[i,3]],InterpolationOrder->0,PlotRangePadding->None,PlotRange->All,PlotLegends->BarLegend[Automatic,LegendLabel->"intensity"],FrameLabel->{"x","y"},PlotLabel->"\[Theta] = "<>ToString@intensityTiltData[[i,1]],ColorFunction->GrayLevel]*)
(*,{i,1,21,1},TrackedSymbols:>{i}]*)


(* ::Subsection:: *)
(*maltese cross as a function of tilt angles \[Theta], \[Phi]*)


(* ::Subsubsection:: *)
(*multiple runs*)


(* ::Input:: *)
(*solution=solveHemisphereLC[0.0,0.0];*)
(*setDefaultParameters[];*)
(*\[Lambda]=4.7`*^-7;*)
(*(*\[Lambda]=5.5`*^-7;*)*)
(*(*\[Lambda]=6.3`*^-7;*)*)
(*setRefractionIndices[\[Lambda]];*)
(*{xVals,yVals,zVals,nxJones,nyJones,nzJones,d}=setJonesCalculusDiscretization[rDropletDim];*)
(**)
(*(*outputDirName="true_compiled_test_dropletTilt_thetaPhi_analyticDirector_wvl_"<>ToString[Round[\[Lambda] 10^9]]<>"_DnLC_"<>ToString[DecimalForm[\[CapitalDelta]nLC,{3,2}]];*)*)
(*outputDirName="gravity_wvl_"<>ToString[Round[\[Lambda] 10^9]]<>"_DnLC_"<>ToString[DecimalForm[\[CapitalDelta]nLC,{3,2}]];*)
(*outputFileFormats={"h5"(*,"mat"*)};*)
(*plotQ=True;*)
(*compiledQ=False;*)
(**)
(*(** parameter scan **)*)
(*(*{nThetas,nPhis}={200,100};*)*)
(*{nThetas,nPhis}={10,20};*)
(*tiltAnglesTheta=N@Subdivide[0,\[Pi]/2,nThetas-1];*)
(*tiltAnglesPhi=N@Subdivide[0,\[Pi]/4,nPhis-1];*)
(**)
(*(** manual specification of parameters **)*)
(*(*{tiltAnglesTheta,tiltAnglesPhi}=N@{(*0.99\[Pi]*){0},{0}};*)*)
(*(*tiltAnglesTheta={0,0.001,0.01,0.1};*)*)
(*(*tiltAnglesTheta=N@{0.1\[Pi]/2};*)*)
(*tiltAnglesPhi=N@{0};*)
(*(*tiltAnglesTheta=N@{\[Pi]/32,\[Pi]/16,\[Pi]/8};*)
(*tiltAnglesPhi=N@{0,\[Pi]/8,\[Pi]/4};*)*)
(*{nThetas,nPhis}=Length/@{tiltAnglesTheta,tiltAnglesPhi};*)
(**)
(*(** progress indicator **)*)
(*progressIndicatorRemainingTime[];*)
(*If[$VersionNumber<13,outputFileFormats=DeleteCases[outputFileFormats,"mat"]];*)
(**)
(*(** create data **)*)
(*timing=AbsoluteTiming[*)
(*intensityTiltData=Table[*)
(*directorAnglesRefIndices3d=calcLCn\[Phi][tiltAngleTheta,tiltAnglePhi,xVals,yVals,zVals,compiledQ];*)
(**)
(*(** calculate polarization state of light passing through water, LC and FC **)*)
(*afterDropletEField=jonesCalculusPropagator[xVals,yVals,directorAnglesRefIndices3d,nLCo,d,rDropletDim,nxJones,nyJones,nzJones,\[Lambda],compiledQ];*)
(**)
(*intensityXY=Map[Abs[#[[1]]]^2+Abs[#[[2]]]^2&,afterDropletEField,{2}];*)
(*avgIntensity=Mean[Flatten[intensityXY]];*)
(**)
(*If[plotQ,*)
(*plot=ListDensityPlot[intensityXY,InterpolationOrder->0,PlotRangePadding->None,PlotRange->All,PlotLegends->BarLegend[Automatic,LegendLabel->"intensity"],FrameLabel->{"x","y"},PlotLabel->"\[Theta] = "<>ToString[DecimalForm[tiltAngleTheta 180/\[Pi],{4,1}]]<>"\[ThinSpace]\[Degree]"*)
(*,ColorFunction->GrayLevel];*)
(*Export[FileNameJoin[{NotebookDirectory[],outputDirName,"afterDropletIntensityXY_"<>ToString[Round[\[Lambda] 10^9]]<>"_theta_"<>ToString[DecimalForm[tiltAngleTheta,{4,3}]]<>"_phi_"<>ToString[DecimalForm[tiltAnglePhi,{4,3}]]<>".png"}],plot];*)
(*];*)
(**)
(*(** save data **)*)
(*finalEfieldState=Map[ReIm,afterDropletEField];*)
(*Export[FileNameJoin[{NotebookDirectory[],outputDirName,"afterDropletEField_"<>ToString[Round[\[Lambda] 10^9]]<>"_theta_"<>ToString[DecimalForm[tiltAngleTheta,{4,3}]]<>"_phi_"<>ToString[DecimalForm[tiltAnglePhi,{4,3}]]<>"."<>#}],finalEfieldState]&/@outputFileFormats;*)
(**)
(*(** update progress indicator **)*)
(*updateProgressIndicator[Length[tiltAnglesTheta]Length[tiltAnglesPhi]];*)
(**)
(*(** output **)*)
(*{tiltAngleTheta,tiltAnglePhi,avgIntensity,afterDropletEField,intensityXY}*)
(*,{tiltAngleTheta,tiltAnglesTheta}*)
(*,{tiltAnglePhi,tiltAnglesPhi}];*)
(*][[1]];*)
(*If[NumericQ[timing],Print[If[timing>3600,ToString[timing/3600]<>"\[ThinSpace]h",IntegerString[Round[timing]]<>"\[ThinSpace]s"]]];*)


(* ::Subsubsection::Closed:: *)
(*check optical pathlength*)


(* ::Input:: *)
(*opticalPathLength=7000/nzJones Table[Total[directorAnglesRefIndices3d[[i,j,All,2]]],{i,nxJones},{j,nyJones}];*)
(*opticalPathLength[[1,1]]*)
(*3500nLCe+3500nOil*)
(*3500nLCo+3500nOil*)
(*7000nWater*)
(*7000nLCe*)
(**)
(*ListDensityPlot[opticalPathLength,PlotLegends->Automatic,PlotRangePadding->None,FrameStyle->Directive[Black,18],FrameLabel->{"x","z"}]*)
(*ListDensityPlot[Mod[opticalPathLength,\[Lambda] 10^9],PlotLegends->Automatic,PlotRangePadding->None,FrameStyle->Directive[Black,18],FrameLabel->{"x","z"}]*)


(* ::Subsubsection::Closed:: *)
(*single run*)


(* ::Input:: *)
(*{tiltAngleTheta,tiltAnglePhi}=N@{\[Pi]/4,\[Pi]/8};*)


(* ::Input:: *)
(*rotaMatrix=RotationMatrix[{{0,0,1},tiltDirector[tiltAngleTheta,tiltAnglePhi]}];*)
(*rotaMatrixInv=RotationMatrix[{tiltDirector[tiltAngleTheta,tiltAnglePhi],{0,0,1}}];*)


(* ::Input:: *)
(*(** LC director projected on xy plane angle = marker for liquid type (LC,FC,water) **)*)
(*(** calculate extraordinary refractive index Subscript[n, e]' for LC director (projected length normal to wave vector k) **)*)
(*directorAnglesRefIndices3d=Table[*)
(*position={x,y,z};*)
(*r=Norm[position];*)
(*{directorAngle,refIndexE}=If[r<=0.99,*)
(*(** inside droplet **)*)
(*If[position . tiltDirector[tiltAngleTheta,tiltAnglePhi]>0.0,*)
(*rotatedPosition=rotaMatrixInv . position;*)
(*rotatedDirector=rotaMatrix . director[rotatedPosition];*)
(*lcAngle=ArcTan[rotatedDirector[[1]],rotatedDirector[[2]]];*)
(*nLCeMod=calcXYProjectedDirectorLength[rotatedDirector,waveVector,nLCe];*)
(*{lcAngle,nLCeMod},{noAngle,nOil}]*)
(*(** outside droplet **)*)
(*,{noAngle,nWater}]*)
(*,{x,xVals},{y,yVals},{z,zVals}];*)


(* ::Input:: *)
(*(** calculate polarization state of light in LC part + FC part **)*)
(*afterDropletEField=Table[*)
(*directorAnglesRefIndices=Reverse@directorAnglesRefIndices3d[[xIndex,yIndex]];*)
(*singleJonesMatrices=Map[( *)
(*{directorAngle,refIndex}=#;*)
(*If[-\[Pi]<=directorAngle<=\[Pi],*)
(*jonesMatrixLC[directorAngle,\[Delta]LC[refIndex-nLCo],\[Delta]Iso[nLCo]],*)
(*jonesMatrixIsotropic[\[Delta]Iso[refIndex]]*)
(*]*)
(*)&,directorAnglesRefIndices];*)
(*dropletJonesMatrix=Dot@@singleJonesMatrices;*)
(*afterDropletState=dropletJonesMatrix . polarizationStateHorizontal;*)
(*jonesMatrixVerticalPolarizer . afterDropletState*)
(*,{xIndex,nxJones}*)
(*,{yIndex,nyJones}*)
(*];*)
(*(*{phaseX,phaseY}=Arg/@afterDropletEField;*)*)
(*intensity=Map[Abs[#[[1]]]^2+Abs[#[[2]]]^2&,afterDropletEField,{2}];*)
(*avgIntensity=Mean[Flatten[intensity]];*)


(* ::Input:: *)
(*(** save data **)*)
(*finalEfieldState=Map[ReIm,afterDropletEField];*)
(*pthout=FileNameJoin[{NotebookDirectory[],"dropletTilt_thetaPhi_lcfch2o_5"}];*)
(*Export[FileNameJoin[{pthout,"afterDropletEField_550_theta_"<>ToString[DecimalForm[tiltAngleTheta,{4,3}]]<>"_phi_"<>ToString[DecimalForm[tiltAnglePhi,{4,3}]]<>".h5"}],finalEfieldState];*)
(**)
(*(** update progress indicator **)*)
(*pi=N[++counter/(Length[tiltAnglesTheta]Length[tiltAnglesPhi])];*)
(**)
(*(** output **)*)
(*output={tiltAngleTheta,tiltAnglePhi,avgIntensity,afterDropletEField,intensity};*)


(* ::Subsubsection::Closed:: *)
(*save E fields*)


(* ::Input:: *)
(*(** save output electric field **)*)
(*saveStateQ=True;*)
(*If[saveStateQ,*)
(*Do[*)
(*finalEfieldState=Map[ReIm,intensityTiltData[[i,j,4]]];*)
(*Export[FileNameJoin[{$HomeDirectory,"Desktop","dropletTiltTest_brightSpot","afterDropletEField_550_i_"<>ToString[i]<>"_j_"<>ToString[j]<>".h5"}],finalEfieldState];*)
(*,{i,Length[tiltAnglesTheta]},{j,Length[tiltAnglesTheta]}]*)
(*];*)


(* ::Subsubsection:: *)
(*show intensity*)


(* ::Input:: *)
(*intensityX=Map[Abs[#[[1]]]^2&,afterDropletEField,{2}];*)
(*intensityY=Map[Abs[#[[2]]]^2&,afterDropletEField,{2}];*)


(* ::Input:: *)
(*i=-1;*)
(*plot=ListDensityPlot[intensityX,InterpolationOrder->0,PlotRangePadding->None,PlotRange->All,PlotLegends->BarLegend[Automatic,LegendLabel->"intensity"],FrameLabel->{"x","y"},PlotLabel->"\[Theta] = "<>ToString[DecimalForm[intensityX 180/\[Pi],{4,1}]]<>" \[Degree]"*)
(*,ColorFunction->GrayLevel*)
(*]*)


(* ::Input:: *)
(*i=-1;*)
(*plot=ListDensityPlot[intensityY,InterpolationOrder->0,PlotRangePadding->None,PlotRange->All,PlotLegends->BarLegend[Automatic,LegendLabel->"intensity"],FrameLabel->{"x","y"},PlotLabel->"\[Theta] = "<>ToString[DecimalForm[intensityXY 180/\[Pi],{4,1}]]<>" \[Degree]"*)
(*,ColorFunction->GrayLevel*)
(*(*,PerformanceGoal->"Speed"*)]*)


(* ::Input:: *)
(*i=-1;*)
(*plot=ListDensityPlot[intensityXY,InterpolationOrder->0,PlotRangePadding->None,PlotRange->All,PlotLegends->BarLegend[Automatic,LegendLabel->"intensity"],FrameLabel->{"x","y"},PlotLabel->"\[Theta] = "<>ToString[DecimalForm[intensityXY 180/\[Pi],{4,1}]]<>" \[Degree]"*)
(*,ColorFunction->GrayLevel*)
(*(*,PerformanceGoal->"Speed"*)]*)


(* ::Input:: *)
(*i=-1;*)
(*plot=ListDensityPlot[intensityTiltData[[i,1,5]],InterpolationOrder->0,PlotRangePadding->None,PlotRange->All,PlotLegends->BarLegend[Automatic,LegendLabel->"intensity"],FrameLabel->{"x","y"},PlotLabel->"\[Theta] = "<>ToString[DecimalForm[intensityTiltData[[i,1,1]]180/\[Pi],{4,1}]]<>" \[Degree]"*)
(*,ColorFunction->GrayLevel*)
(*(*,PerformanceGoal->"Speed"*)]*)


(* ::Input:: *)
(*Manipulate[*)
(*intensityXY=intensityTiltData[[i,j,5]];*)
(*(*ListDensityPlot[intensityTiltData\[LeftDoubleBracket]i,j,5\[RightDoubleBracket],InterpolationOrder->0,PlotRangePadding->None,PlotRange->All,PlotLegends->BarLegend[Automatic,LegendLabel\[Rule]"intensity"],FrameLabel->{"x","y"},PlotLabel->"\[Theta] = "<>ToString[DecimalForm[intensityTiltData\[LeftDoubleBracket]i,j,1\[RightDoubleBracket]180/\[Pi],{4,1}]]<>"\[ThinSpace]\[Degree], "<>"\[Phi] = "<>ToString[DecimalForm[intensityTiltData\[LeftDoubleBracket]i,j,2\[RightDoubleBracket]180/\[Pi],{4,1}]]<>"\[ThinSpace]\[Degree]"*)
(*,ColorFunction->GrayLevel,PerformanceGoal->"Quality"*)
(*(*,PerformanceGoal->"Speed"*)]*)*)
(*ListDensityPlot[intensityXY,InterpolationOrder->0,PlotRangePadding->None,PlotRange->All,PlotLegends->BarLegend[Automatic,LegendLabel->"intensity \!\(\**)
(*StyleBox[SubscriptBox[\"I\", \"xy\"],\nFontSlant->\"Italic\"]\)"],FrameLabel->{"x (\[Mu]m)","y (\[Mu]m)"},FrameStyle->Directive[Black],DataRange->{{-3.5,3.5},{-3.5,3.5}},PlotLabel->Style["\[Theta] = "<>ToString[DecimalForm[intensityTiltData[[i,j,1]]180/\[Pi],{4,1}]]<>"\[ThinSpace]\[Degree], "<>"\[Phi] = "<>ToString[DecimalForm[intensityTiltData[[i,j,2]]180/\[Pi],{4,1}]]<>"\[ThinSpace]\[Degree]",Black]*)
(*,ColorFunction->GrayLevel,PerformanceGoal->"Quality"*)
(*(*,PerformanceGoal->"Speed"*)]*)
(*,{i,1,nThetas,1,Appearance->"Open"}*)
(*,{j,1,nPhis,1,Appearance->"Open"}*)
(*,TrackedSymbols:>{i,j}*)
(*]*)


(* ::Subsubsection:: *)
(*scalar intensities right after droplet - lookup map*)


(* ::Input:: *)
(*plot=ListDensityPlot[{180/\[Pi]#[[1]],180/\[Pi]#[[2]],#[[3]]}&/@Flatten[intensityTiltData[[All,All,{1,2,3}]],1],FrameLabel->{"\[Theta] (\[Degree])","\[Phi] (\[Degree])"},AxesStyle->Black,Mesh->None,PlotRangePadding->None,FrameStyle->Directive[Black],AspectRatio->0.5]*)
(*Export[FileNameJoin[{NotebookDirectory[],outputDirName,"thetaPhi_intensity_2d.png"}],plot]*)


(* ::Input:: *)
(*plot=Show[*)
(*ListPlot3D[{180/\[Pi]#[[1]],180/\[Pi]#[[2]],#[[3]]}&/@Flatten[intensityTiltData[[All,All,{1,2,3}]],1],AxesLabel->{"\[Theta] (\[Degree])","\[Phi] (\[Degree])","I"},AxesStyle->Black,Mesh->None,BoxRatios->{1,0.5,0.5}],*)
(*ListPointPlot3D[{180/\[Pi]#[[1]],180/\[Pi]#[[2]],#[[3]]}&/@Flatten[intensityTiltData[[All,All,{1,2,3}]],1],AxesLabel->{"\[Theta] (\[Degree])","\[Phi] (\[Degree])","I"},AxesStyle->Black,PlotStyle->Black]*)
(*]*)
(*Export[FileNameJoin[{NotebookDirectory[],outputDirName,"thetaPhi_intensity_3d.png"}],plot]*)


(* ::Subsubsection::Closed:: *)
(*create movie of tilted intensities for different resolutions*)


(* ::Input:: *)
(*Manipulate[*)
(*ListDensityPlot[intensityTiltData[[i,j,5]],InterpolationOrder->0,PlotRangePadding->None,PlotRange->All,PlotLegends->BarLegend[Automatic,LegendLabel->"intensity"],FrameLabel->{"x","y"},PlotLabel->"\[Theta] = "<>ToString@intensityTiltData[[i,j,1]],ColorFunction->GrayLevel]*)
(*,{i,1,Length[tiltAnglesTheta],1}*)
(*,{j,1,Length[tiltAnglesPhi],1}*)
(*]*)


(* ::Input:: *)
(*counter=0;*)
(*Do[*)
(*label="\[Theta] = "<>ToString[DecimalForm[intensityTiltData[[i,1]]180/\[Pi],{4,1}]]<>"\[Degree], \[Phi] = "<>ToString[DecimalForm[tiltAnglePhi 180/\[Pi],{4,1}]]<>"\[Degree]";*)
(*plotHighRes=ListDensityPlot[intensityTiltData[[i,j,4]],InterpolationOrder->0,PlotRangePadding->None,PlotRange->All,PlotLegends->BarLegend[Automatic,LegendLabel->"intensity"],FrameLabel->{"x","y"},ColorFunction->GrayLevel,ImageSize->{Automatic,400}];*)
(*plotLowRes=ListDensityPlot[intensityTiltData[[i,j,4]],InterpolationOrder->0,PlotRangePadding->None,PlotRange->All,PlotLegends->BarLegend[Automatic,LegendLabel->"intensity"],FrameLabel->{"x","y"},ColorFunction->GrayLevel,ImageSize->{Automatic,400},PerformanceGoal->"Speed"];*)
(*plot=Grid[{{label,SpanFromLeft},{plotHighRes,plotLowRes}}];*)
(*Export[FileNameJoin[{NotebookDirectory[],"mma_polarizationMovie_"<>ToString[tiltAnglePhi 180/\[Pi]]<>"_deg","p_"<>IntegerString[++counter,10,IntegerLength[Length@intensityTiltData]+1]<>".png"}],plot];*)
(*,{i,Length[intensityTiltData]}];*)


(* ::Subsubsection::Closed:: *)
(*debug*)


(* ::Input:: *)
(*avgRefIndexE=Table[*)
(*Mean[directorAnglesRefIndices3d[[xIndex,yIndex,All,2]]]*)
(*,{xIndex,nx},{yIndex,ny}*)
(*];*)


(* ::Input:: *)
(*ListDensityPlot[avgRefIndexE,InterpolationOrder->0,PlotRange->All,PlotLegends->Automatic]*)


(* ::Input:: *)
(*ListPlot[avgRefIndexE[[50]]]*)


(* ::Input:: *)
(*totalRefIndexE=Table[*)
(*Mod[Total[directorAnglesRefIndices3d[[xIndex,yIndex,All,2]]]dz rDropletDim 10^6,\[Lambda] 10^6]*)
(*,{xIndex,nx},{yIndex,ny}*)
(*];*)
(*ListDensityPlot[totalRefIndexE,InterpolationOrder->0,PlotRange->All,PlotLegends->Automatic]*)


(* ::Input:: *)
(*Histogram[Flatten[directorAnglesRefIndices3d[[All,All,All,3]]]]*)


(* ::Input:: *)
(*\[CapitalDelta]n=directorAnglesRefIndices3d[[50,All,All,2]]\[Transpose]-nLCo;*)


(* ::Input:: *)
(*ListDensityPlot[\[CapitalDelta]n,PlotLegends->Automatic,PlotRange->All,InterpolationOrder->0]*)


(* ::Input:: *)
(*avgnLC=Table[*)
(*Mean[directorAnglesRefIndices3d[[xIndex,yIndex,All,2]]]*)
(*,{xIndex,Dimensions[directorAnglesRefIndices3d][[1]]}*)
(*,{yIndex,Dimensions[directorAnglesRefIndices3d][[2]]}];*)
(*ListDensityPlot[avgnLC,InterpolationOrder->0]*)


(* ::Input:: *)
(*ListContourPlot[avgnLC,Contours->{nz nWater,1400,1410(*nz(0.5 nLCo+0.5nLCe)*)}]*)


(* ::Subsubsection::Closed:: *)
(*old (without water)*)


(* ::Input:: *)
(*(*afterDropletState=dropletJonesMatrix.RotationMatrix[polarizerAnglePhi].polarizationStateHorizontal;*)
(*Inverse[RotationMatrix[polarizerAnglePhi]].jonesMatrixVerticalPolarizer.RotationMatrix[polarizerAnglePhi].afterDropletState*)
(**)
(*directorAnglesRefIndices=Reverse@directorAnglesRefIndices3d[[xIndex,yIndex]];*)
(*filteredAnglesRefIndices=If[#\[LeftDoubleBracket]1\[RightDoubleBracket]!=waterAngle,#,Nothing]&/@directorAnglesRefIndices;*)
(*If[filteredAnglesRefIndices!={},*)
(*(** ray intersection with the droplet **)*)
(*(*singleJonesMatrices=Map[If[#\[LeftDoubleBracket]1\[RightDoubleBracket]!=fcAngle,jonesMatrixLC[#\[LeftDoubleBracket]1\[RightDoubleBracket],nLCe],jonesMatrixFC]&,filteredAnglesRefIndices];*)*)
(*singleJonesMatrices=Map[If[#\[LeftDoubleBracket]1\[RightDoubleBracket]!=fcAngle,jonesMatrixLC[#\[LeftDoubleBracket]1\[RightDoubleBracket],(*nLCe*)#\[LeftDoubleBracket]2\[RightDoubleBracket]],jonesMatrixFC]&,filteredAnglesRefIndices];dropletJonesMatrix=Dot@@singleJonesMatrices;*)
(*afterDropletState=dropletJonesMatrix.polarizationStateHorizontal;*)
(*jonesMatrixVerticalPolarizer.afterDropletState,*)
(*(** no ray intersection with the droplet **)*)
(*jonesMatrixVerticalPolarizer.polarizationStateHorizontal*)
(*]*)*)


(* ::Subsection::Closed:: *)
(*maltese cross as a function of defect position angles \[Theta], \[Phi]*)


(* ::Subsubsection:: *)
(*multiple runs*)


(* ::Input:: *)
(*(*\[CapitalDelta]nLC=0.83;*)
(*setRefractionIndices[\[Lambda],\[CapitalDelta]nLC];*)*)
(*setDefaultParameters[];*)
(*\[Lambda]=4.7`*^-7;*)
(*(*\[Lambda]=5.5`*^-7;*)*)
(*(*\[Lambda]=6.3`*^-7;*)*)
(*setRefractionIndices[\[Lambda]];*)
(*{xVals,yVals,zVals,nxJones,nyJones,nzJones,d}=setJonesCalculusDiscretization[rDropletDim];*)
(**)
(*outputDirName="displaced_defect_grav_tilt_wvl_"<>ToString[Round[\[Lambda] 10^9]]<>"_DnLC_"<>ToString[DecimalForm[\[CapitalDelta]nLC,{3,2}]];*)
(*outputFileFormats={"h5"(*,"mat"*)};*)
(*plotQ=True;*)
(*compiledQ=False;*)
(**)
(*(** parameter scan **)*)
(*(*{nThetas,nPhis}={200,100};*)*)
(*{nThetas,nPhis}={40,20};*)
(*(*nThetas=40;*)*)
(*tiltAnglesTheta=N@Subdivide[0,\[Pi]/2,nThetas-1];*)
(*tiltAnglesPhi=N@Subdivide[0,\[Pi]/4,nPhis-1];*)
(**)
(*(** manual specification of parameters **)*)
(*(*{tiltAnglesTheta,tiltAnglesPhi}=N@{(*0.99\[Pi]*){0},{0}};*)*)
(*(*tiltAnglesTheta={0,0.001,0.01,0.1};*)*)
(*(*tiltAnglesTheta=N@{0.1\[Pi]/2};*)*)
(*tiltAnglesPhi=N@{3\[Pi]/2};*)
(*(*tiltAnglesTheta=N@{\[Pi]/32,\[Pi]/16,\[Pi]/8};*)
(*tiltAnglesPhi=N@{0,\[Pi]/8,\[Pi]/4};*)*)
(*{nThetas,nPhis}=Length/@{tiltAnglesTheta,tiltAnglesPhi};*)
(**)
(*(** progress indicator **)*)
(*progressIndicatorRemainingTime[];*)
(*If[$VersionNumber<13,outputFileFormats=DeleteCases[outputFileFormats,"mat"]];*)
(**)
(*(** create data **)*)
(*timing=AbsoluteTiming[*)
(*intensityTiltData=Table[*)
(*solution=solveHemisphereLC[0.5tiltAngleTheta,0.0];*)
(**)
(*directorAnglesRefIndices3d=calcLCn\[Phi][0.5tiltAngleTheta,tiltAnglePhi,xVals,yVals,zVals,compiledQ];*)
(**)
(*(** calculate polarization state of light passing through water, LC and FC **)*)
(*afterDropletEField=jonesCalculusPropagator[xVals,yVals,directorAnglesRefIndices3d,nLCo,d,rDropletDim,nxJones,nyJones,nzJones,\[Lambda],compiledQ];*)
(**)
(*intensityXY=Map[Abs[#[[1]]]^2+Abs[#[[2]]]^2&,afterDropletEField,{2}];*)
(*avgIntensity=Mean[Flatten[intensityXY]];*)
(**)
(*If[plotQ,*)
(*plot=ListDensityPlot[intensityXY,InterpolationOrder->0,PlotRangePadding->None,PlotRange->All,PlotLegends->BarLegend[Automatic,LegendLabel->"intensity"],FrameLabel->{"x","y"},PlotLabel->"\[Theta] = "<>ToString[DecimalForm[tiltAngleTheta 180/\[Pi],{4,1}]]<>"\[ThinSpace]\[Degree]"*)
(*,ColorFunction->GrayLevel];*)
(*Export[FileNameJoin[{NotebookDirectory[],outputDirName,"afterDropletIntensityXY_"<>ToString[Round[\[Lambda] 10^9]]<>"_theta_"<>ToString[DecimalForm[tiltAngleTheta,{4,3}]]<>"_phi_"<>ToString[DecimalForm[tiltAnglePhi,{4,3}]]<>".png"}],plot];*)
(*];*)
(**)
(*(** save data **)*)
(*finalEfieldState=Map[ReIm,afterDropletEField];*)
(*Export[FileNameJoin[{NotebookDirectory[],outputDirName,"afterDropletEField_"<>ToString[Round[\[Lambda] 10^9]]<>"_theta_"<>ToString[DecimalForm[tiltAngleTheta,{4,3}]]<>"_phi_"<>ToString[DecimalForm[tiltAnglePhi,{4,3}]]<>"."<>#}],finalEfieldState]&/@outputFileFormats;*)
(**)
(*(** update progress indicator **)*)
(*updateProgressIndicator[Length[tiltAnglesTheta]Length[tiltAnglesPhi]];*)
(**)
(*(** output **)*)
(*{tiltAngleTheta,tiltAnglePhi,avgIntensity,afterDropletEField,intensityXY}*)
(*,{tiltAngleTheta,tiltAnglesTheta}*)
(*,{tiltAnglePhi,tiltAnglesPhi}*)
(*];*)
(*][[1]];*)
(*If[NumericQ[timing],Print[If[timing>3600,ToString[timing/3600]<>"\[ThinSpace]h",IntegerString[Round[timing]]<>"\[ThinSpace]s"]]];*)


(* ::Subsubsection::Closed:: *)
(*check optical pathlength*)


(* ::Input:: *)
(*opticalPathLength=7000/nzJones Table[Total[directorAnglesRefIndices3d[[i,j,All,2]]],{i,nxJones},{j,nyJones}];*)
(*opticalPathLength[[1,1]]*)
(*3500nLCe+3500nOil*)
(*3500nLCo+3500nOil*)
(*7000nWater*)
(*7000nLCe*)
(**)
(*ListDensityPlot[opticalPathLength,PlotLegends->Automatic,PlotRangePadding->None,FrameStyle->Directive[Black,18],FrameLabel->{"x","z"}]*)
(*ListDensityPlot[Mod[opticalPathLength,\[Lambda] 10^9],PlotLegends->Automatic,PlotRangePadding->None,FrameStyle->Directive[Black,18],FrameLabel->{"x","z"}]*)


(* ::Subsubsection::Closed:: *)
(*single run*)


(* ::Input:: *)
(*{tiltAngleTheta,tiltAnglePhi}=N@{\[Pi]/4,\[Pi]/8};*)


(* ::Input:: *)
(*rotaMatrix=RotationMatrix[{{0,0,1},tiltDirector[tiltAngleTheta,tiltAnglePhi]}];*)
(*rotaMatrixInv=RotationMatrix[{tiltDirector[tiltAngleTheta,tiltAnglePhi],{0,0,1}}];*)


(* ::Input:: *)
(*(** LC director projected on xy plane angle = marker for liquid type (LC,FC,water) **)*)
(*(** calculate extraordinary refractive index Subscript[n, e]' for LC director (projected length normal to wave vector k) **)*)
(*directorAnglesRefIndices3d=Table[*)
(*position={x,y,z};*)
(*r=Norm[position];*)
(*{directorAngle,refIndexE}=If[r<=0.99,*)
(*(** inside droplet **)*)
(*If[position . tiltDirector[tiltAngleTheta,tiltAnglePhi]>0.0,*)
(*rotatedPosition=rotaMatrixInv . position;*)
(*rotatedDirector=rotaMatrix . director[rotatedPosition];*)
(*lcAngle=ArcTan[rotatedDirector[[1]],rotatedDirector[[2]]];*)
(*nLCeMod=calcXYProjectedDirectorLength[rotatedDirector,waveVector,nLCe];*)
(*{lcAngle,nLCeMod},{noAngle,nOil}]*)
(*(** outside droplet **)*)
(*,{noAngle,nWater}]*)
(*,{x,xVals},{y,yVals},{z,zVals}];*)


(* ::Input:: *)
(*(** calculate polarization state of light in LC part + FC part **)*)
(*afterDropletEField=Table[*)
(*directorAnglesRefIndices=Reverse@directorAnglesRefIndices3d[[xIndex,yIndex]];*)
(*singleJonesMatrices=Map[( *)
(*{directorAngle,refIndex}=#;*)
(*If[-\[Pi]<=directorAngle<=\[Pi],*)
(*jonesMatrixLC[directorAngle,\[Delta]LC[refIndex-nLCo],\[Delta]Iso[nLCo]],*)
(*jonesMatrixIsotropic[\[Delta]Iso[refIndex]]*)
(*]*)
(*)&,directorAnglesRefIndices];*)
(*dropletJonesMatrix=Dot@@singleJonesMatrices;*)
(*afterDropletState=dropletJonesMatrix . polarizationStateHorizontal;*)
(*jonesMatrixVerticalPolarizer . afterDropletState*)
(*,{xIndex,nxJones}*)
(*,{yIndex,nyJones}*)
(*];*)
(*(*{phaseX,phaseY}=Arg/@afterDropletEField;*)*)
(*intensity=Map[Abs[#[[1]]]^2+Abs[#[[2]]]^2&,afterDropletEField,{2}];*)
(*avgIntensity=Mean[Flatten[intensity]];*)


(* ::Input:: *)
(*(** save data **)*)
(*finalEfieldState=Map[ReIm,afterDropletEField];*)
(*pthout=FileNameJoin[{NotebookDirectory[],"dropletTilt_thetaPhi_lcfch2o_5"}];*)
(*Export[FileNameJoin[{pthout,"afterDropletEField_550_theta_"<>ToString[DecimalForm[tiltAngleTheta,{4,3}]]<>"_phi_"<>ToString[DecimalForm[tiltAnglePhi,{4,3}]]<>".h5"}],finalEfieldState];*)
(**)
(*(** update progress indicator **)*)
(*pi=N[++counter/(Length[tiltAnglesTheta]Length[tiltAnglesPhi])];*)
(**)
(*(** output **)*)
(*output={tiltAngleTheta,tiltAnglePhi,avgIntensity,afterDropletEField,intensity};*)


(* ::Subsubsection::Closed:: *)
(*save E fields*)


(* ::Input:: *)
(*(** save output electric field **)*)
(*saveStateQ=True;*)
(*If[saveStateQ,*)
(*Do[*)
(*finalEfieldState=Map[ReIm,intensityTiltData[[i,j,4]]];*)
(*Export[FileNameJoin[{$HomeDirectory,"Desktop","dropletTiltTest_brightSpot","afterDropletEField_550_i_"<>ToString[i]<>"_j_"<>ToString[j]<>".h5"}],finalEfieldState];*)
(*,{i,Length[tiltAnglesTheta]},{j,Length[tiltAnglesTheta]}]*)
(*];*)


(* ::Subsubsection:: *)
(*show intensity*)


(* ::Input:: *)
(*intensityX=Map[Abs[#[[1]]]^2&,afterDropletEField,{2}];*)
(*intensityY=Map[Abs[#[[2]]]^2&,afterDropletEField,{2}];*)


(* ::Input:: *)
(*i=-1;*)
(*plot=ListDensityPlot[intensityX,InterpolationOrder->0,PlotRangePadding->None,PlotRange->All,PlotLegends->BarLegend[Automatic,LegendLabel->"intensity"],FrameLabel->{"x","y"},PlotLabel->"\[Theta] = "<>ToString[DecimalForm[intensityX 180/\[Pi],{4,1}]]<>" \[Degree]"*)
(*,ColorFunction->GrayLevel*)
(*]*)


(* ::Input:: *)
(*i=-1;*)
(*plot=ListDensityPlot[intensityY,InterpolationOrder->0,PlotRangePadding->None,PlotRange->All,PlotLegends->BarLegend[Automatic,LegendLabel->"intensity"],FrameLabel->{"x","y"},PlotLabel->"\[Theta] = "<>ToString[DecimalForm[intensityXY 180/\[Pi],{4,1}]]<>" \[Degree]"*)
(*,ColorFunction->GrayLevel*)
(*(*,PerformanceGoal->"Speed"*)]*)


(* ::Input:: *)
(*i=-1;*)
(*plot=ListDensityPlot[intensityXY,InterpolationOrder->0,PlotRangePadding->None,PlotRange->All,PlotLegends->BarLegend[Automatic,LegendLabel->"intensity"],FrameLabel->{"x","y"},PlotLabel->"\[Theta] = "<>ToString[DecimalForm[intensityXY 180/\[Pi],{4,1}]]<>" \[Degree]"*)
(*,ColorFunction->GrayLevel*)
(*(*,PerformanceGoal->"Speed"*)]*)


(* ::Input:: *)
(*i=-1;*)
(*plot=ListDensityPlot[intensityTiltData[[i,1,5]],InterpolationOrder->0,PlotRangePadding->None,PlotRange->All,PlotLegends->BarLegend[Automatic,LegendLabel->"intensity"],FrameLabel->{"x","y"},PlotLabel->"\[Theta] = "<>ToString[DecimalForm[intensityTiltData[[i,1,1]]180/\[Pi],{4,1}]]<>" \[Degree]"*)
(*,ColorFunction->GrayLevel*)
(*(*,PerformanceGoal->"Speed"*)]*)


(* ::Input:: *)
(*Manipulate[*)
(*intensityXY=intensityTiltData[[i,j,5]];*)
(*(*ListDensityPlot[intensityTiltData\[LeftDoubleBracket]i,j,5\[RightDoubleBracket],InterpolationOrder->0,PlotRangePadding->None,PlotRange->All,PlotLegends->BarLegend[Automatic,LegendLabel\[Rule]"intensity"],FrameLabel->{"x","y"},PlotLabel->"\[Theta] = "<>ToString[DecimalForm[intensityTiltData\[LeftDoubleBracket]i,j,1\[RightDoubleBracket]180/\[Pi],{4,1}]]<>"\[ThinSpace]\[Degree], "<>"\[Phi] = "<>ToString[DecimalForm[intensityTiltData\[LeftDoubleBracket]i,j,2\[RightDoubleBracket]180/\[Pi],{4,1}]]<>"\[ThinSpace]\[Degree]"*)
(*,ColorFunction->GrayLevel,PerformanceGoal->"Quality"*)
(*(*,PerformanceGoal->"Speed"*)]*)*)
(*ListDensityPlot[intensityXY,InterpolationOrder->0,PlotRangePadding->None,PlotRange->All,PlotLegends->BarLegend[Automatic,LegendLabel->"intensity \!\(\**)
(*StyleBox[SubscriptBox[\"I\", \"xy\"],\nFontSlant->\"Italic\"]\)"],FrameLabel->{"x (\[Mu]m)","y (\[Mu]m)"},FrameStyle->Directive[Black],DataRange->{{-3.5,3.5},{-3.5,3.5}},PlotLabel->Style["\[Theta] = "<>ToString[DecimalForm[intensityTiltData[[i,j,1]]180/\[Pi],{4,1}]]<>"\[ThinSpace]\[Degree], "<>"\[Phi] = "<>ToString[DecimalForm[intensityTiltData[[i,j,2]]180/\[Pi],{4,1}]]<>"\[ThinSpace]\[Degree]",Black]*)
(*,ColorFunction->GrayLevel,PerformanceGoal->"Quality"*)
(*(*,PerformanceGoal->"Speed"*)]*)
(*,{i,1,nThetas,1,Appearance->"Open"}*)
(*,{j,1,nPhis,1,Appearance->"Open"}*)
(*,TrackedSymbols:>{i,j}*)
(*]*)


(* ::Subsubsection:: *)
(*scalar intensities right after droplet - lookup map*)


(* ::Input:: *)
(*plot=ListDensityPlot[{180/\[Pi]#[[1]],180/\[Pi]#[[2]],#[[3]]}&/@Flatten[intensityTiltData[[All,All,{1,2,3}]],1],FrameLabel->{"\[Theta] (\[Degree])","\[Phi] (\[Degree])"},AxesStyle->Black,Mesh->None,PlotRangePadding->None,FrameStyle->Directive[Black],AspectRatio->0.5]*)
(*Export[FileNameJoin[{NotebookDirectory[],outputDirName,"thetaPhi_intensity_2d.png"}],plot]*)


(* ::Input:: *)
(*plot=Show[*)
(*ListPlot3D[{180/\[Pi]#[[1]],180/\[Pi]#[[2]],#[[3]]}&/@Flatten[intensityTiltData[[All,All,{1,2,3}]],1],AxesLabel->{"\[Theta] (\[Degree])","\[Phi] (\[Degree])","I"},AxesStyle->Black,Mesh->None,BoxRatios->{1,0.5,0.5}],*)
(*ListPointPlot3D[{180/\[Pi]#[[1]],180/\[Pi]#[[2]],#[[3]]}&/@Flatten[intensityTiltData[[All,All,{1,2,3}]],1],AxesLabel->{"\[Theta] (\[Degree])","\[Phi] (\[Degree])","I"},AxesStyle->Black,PlotStyle->Black]*)
(*]*)
(*Export[FileNameJoin[{NotebookDirectory[],outputDirName,"thetaPhi_intensity_3d.png"}],plot]*)


(* ::Subsubsection::Closed:: *)
(*create movie of tilted intensities for different resolutions*)


(* ::Input:: *)
(*Manipulate[*)
(*ListDensityPlot[intensityTiltData[[i,j,5]],InterpolationOrder->0,PlotRangePadding->None,PlotRange->All,PlotLegends->BarLegend[Automatic,LegendLabel->"intensity"],FrameLabel->{"x","y"},PlotLabel->"\[Theta] = "<>ToString@intensityTiltData[[i,j,1]],ColorFunction->GrayLevel]*)
(*,{i,1,Length[tiltAnglesTheta],1}*)
(*,{j,1,Length[tiltAnglesPhi],1}*)
(*]*)


(* ::Input:: *)
(*counter=0;*)
(*Do[*)
(*label="\[Theta] = "<>ToString[DecimalForm[intensityTiltData[[i,1]]180/\[Pi],{4,1}]]<>"\[Degree], \[Phi] = "<>ToString[DecimalForm[tiltAnglePhi 180/\[Pi],{4,1}]]<>"\[Degree]";*)
(*plotHighRes=ListDensityPlot[intensityTiltData[[i,j,4]],InterpolationOrder->0,PlotRangePadding->None,PlotRange->All,PlotLegends->BarLegend[Automatic,LegendLabel->"intensity"],FrameLabel->{"x","y"},ColorFunction->GrayLevel,ImageSize->{Automatic,400}];*)
(*plotLowRes=ListDensityPlot[intensityTiltData[[i,j,4]],InterpolationOrder->0,PlotRangePadding->None,PlotRange->All,PlotLegends->BarLegend[Automatic,LegendLabel->"intensity"],FrameLabel->{"x","y"},ColorFunction->GrayLevel,ImageSize->{Automatic,400},PerformanceGoal->"Speed"];*)
(*plot=Grid[{{label,SpanFromLeft},{plotHighRes,plotLowRes}}];*)
(*Export[FileNameJoin[{NotebookDirectory[],"mma_polarizationMovie_"<>ToString[tiltAnglePhi 180/\[Pi]]<>"_deg","p_"<>IntegerString[++counter,10,IntegerLength[Length@intensityTiltData]+1]<>".png"}],plot];*)
(*,{i,Length[intensityTiltData]}];*)


(* ::Subsubsection::Closed:: *)
(*debug*)


(* ::Input:: *)
(*avgRefIndexE=Table[*)
(*Mean[directorAnglesRefIndices3d[[xIndex,yIndex,All,2]]]*)
(*,{xIndex,nx},{yIndex,ny}*)
(*];*)


(* ::Input:: *)
(*ListDensityPlot[avgRefIndexE,InterpolationOrder->0,PlotRange->All,PlotLegends->Automatic]*)


(* ::Input:: *)
(*ListPlot[avgRefIndexE[[50]]]*)


(* ::Input:: *)
(*totalRefIndexE=Table[*)
(*Mod[Total[directorAnglesRefIndices3d[[xIndex,yIndex,All,2]]]dz rDropletDim 10^6,\[Lambda] 10^6]*)
(*,{xIndex,nx},{yIndex,ny}*)
(*];*)
(*ListDensityPlot[totalRefIndexE,InterpolationOrder->0,PlotRange->All,PlotLegends->Automatic]*)


(* ::Input:: *)
(*Histogram[Flatten[directorAnglesRefIndices3d[[All,All,All,3]]]]*)


(* ::Input:: *)
(*\[CapitalDelta]n=directorAnglesRefIndices3d[[50,All,All,2]]\[Transpose]-nLCo;*)


(* ::Input:: *)
(*ListDensityPlot[\[CapitalDelta]n,PlotLegends->Automatic,PlotRange->All,InterpolationOrder->0]*)


(* ::Input:: *)
(*avgnLC=Table[*)
(*Mean[directorAnglesRefIndices3d[[xIndex,yIndex,All,2]]]*)
(*,{xIndex,Dimensions[directorAnglesRefIndices3d][[1]]}*)
(*,{yIndex,Dimensions[directorAnglesRefIndices3d][[2]]}];*)
(*ListDensityPlot[avgnLC,InterpolationOrder->0]*)


(* ::Input:: *)
(*ListContourPlot[avgnLC,Contours->{nz nWater,1400,1410(*nz(0.5 nLCo+0.5nLCe)*)}]*)


(* ::Subsubsection::Closed:: *)
(*old (without water)*)


(* ::Input:: *)
(*(*afterDropletState=dropletJonesMatrix.RotationMatrix[polarizerAnglePhi].polarizationStateHorizontal;*)
(*Inverse[RotationMatrix[polarizerAnglePhi]].jonesMatrixVerticalPolarizer.RotationMatrix[polarizerAnglePhi].afterDropletState*)
(**)
(*directorAnglesRefIndices=Reverse@directorAnglesRefIndices3d[[xIndex,yIndex]];*)
(*filteredAnglesRefIndices=If[#\[LeftDoubleBracket]1\[RightDoubleBracket]!=waterAngle,#,Nothing]&/@directorAnglesRefIndices;*)
(*If[filteredAnglesRefIndices!={},*)
(*(** ray intersection with the droplet **)*)
(*(*singleJonesMatrices=Map[If[#\[LeftDoubleBracket]1\[RightDoubleBracket]!=fcAngle,jonesMatrixLC[#\[LeftDoubleBracket]1\[RightDoubleBracket],nLCe],jonesMatrixFC]&,filteredAnglesRefIndices];*)*)
(*singleJonesMatrices=Map[If[#\[LeftDoubleBracket]1\[RightDoubleBracket]!=fcAngle,jonesMatrixLC[#\[LeftDoubleBracket]1\[RightDoubleBracket],(*nLCe*)#\[LeftDoubleBracket]2\[RightDoubleBracket]],jonesMatrixFC]&,filteredAnglesRefIndices];dropletJonesMatrix=Dot@@singleJonesMatrices;*)
(*afterDropletState=dropletJonesMatrix.polarizationStateHorizontal;*)
(*jonesMatrixVerticalPolarizer.afterDropletState,*)
(*(** no ray intersection with the droplet **)*)
(*jonesMatrixVerticalPolarizer.polarizationStateHorizontal*)
(*]*)*)


(* ::Subsection::Closed:: *)
(*maltese cross as a function of tilt angles \[Theta], \[Phi], \[Lambda]*)


(* ::Subsubsection:: *)
(*multiple runs (no mma internal output)*)


(* ::Input:: *)
(*setDefaultParameters[];*)
(*outputFileFormats={"h5"};*)
(*plotQ=False;*)
(**)
(*(** parameter scan **)*)
(*wavelengths=N@Subdivide[470,700,10-1];*)
(*{nThetas,nPhis}={40,20};*)
(*tiltAnglesTheta=N@Subdivide[0,\[Pi]/2,nThetas-1];*)
(*tiltAnglesPhi=N@Subdivide[0,\[Pi]/4,nPhis-1];*)
(**)
(*(** progress indicator **)*)
(*progressIndicatorRemainingTime[];*)
(*If[$VersionNumber<13,outputFileFormats=DeleteCases[outputFileFormats,"mat"]];*)
(**)
(*(** create data **)*)
(*timing=AbsoluteTiming[*)
(*Do[*)
(*\[Lambda]=\[Lambda]nm*10^-9;  (** in m **)*)
(*setRefractionIndices[\[Lambda]];*)
(*{xVals,yVals,zVals,nxJones,nyJones,nzJones,d}=setJonesCalculusDiscretization[rDropletDim];*)
(*outputDirName="compiled_test_dropletTilt_thetaPhi_analyticDirector_wvl_"<>ToString[Round[\[Lambda] 10^9]]<>"_DnLC_"<>ToString[DecimalForm[\[CapitalDelta]nLC,{3,2}]];*)
(*(*psfMicroscopeIncoherent=jinc[xMat,yMat,2na/(\[Lambda] magnification),0,0]^2;*)*)
(**)
(*Do[*)
(*directorAnglesRefIndices3d=calcLCn\[Phi][tiltAngleTheta,tiltAnglePhi,xVals,yVals,zVals,True];*)
(**)
(*(** calculate polarization state of light passing through water, LC and FC **)*)
(*afterDropletEField=jonesCalculusPropagator[xVals,yVals,directorAnglesRefIndices3d,nLCo,d,rDropletDim,nxJones,nyJones,nzJones,\[Lambda]];*)
(**)
(*intensityXY=Map[Abs[#[[1]]]^2+Abs[#[[2]]]^2&,afterDropletEField,{2}];*)
(*intensityXYMicroscopeOptics=ListConvolve[intensityXY,PSFMicroscopeIncoherent];*)
(*avgIntensity=Mean[Flatten[intensityXYMicroscopeOptics]];*)
(**)
(*If[plotQ,*)
(*plot=ListDensityPlot[intensityXY,InterpolationOrder->0,PlotRangePadding->None,PlotRange->All,PlotLegends->BarLegend[Automatic,LegendLabel->"intensity"],FrameLabel->{"x","y"},PlotLabel->"\[Theta] = "<>ToString[DecimalForm[tiltAngleTheta 180/\[Pi],{4,1}]]<>"\[ThinSpace]\[Degree]"*)
(*,ColorFunction->GrayLevel];*)
(*Export[FileNameJoin[{NotebookDirectory[],outputDirName,"afterDropletIntensityXY_"<>ToString[Round[\[Lambda] 10^9]]<>"_theta_"<>ToString[DecimalForm[tiltAngleTheta,{4,3}]]<>"_phi_"<>ToString[DecimalForm[tiltAnglePhi,{4,3}]]<>".png"}],plot];*)
(*];*)
(**)
(*(** save data **)*)
(*finalEfieldState=Map[ReIm,afterDropletEField];*)
(*Export[FileNameJoin[{NotebookDirectory[],outputDirName,"afterDropletEField_"<>ToString[Round[\[Lambda] 10^9]]<>"_theta_"<>ToString[DecimalForm[tiltAngleTheta,{4,3}]]<>"_phi_"<>ToString[DecimalForm[tiltAnglePhi,{4,3}]]<>"."<>#}],finalEfieldState]&/@outputFileFormats;*)
(**)
(*,{tiltAngleTheta,tiltAnglesTheta}*)
(*,{tiltAnglePhi,tiltAnglesPhi}];*)
(**)
(*(** update progress indicator **)*)
(*updateProgressIndicator[Length[tiltAnglesTheta]Length[tiltAnglesPhi]Length[wavelengths]];*)
(*,{\[Lambda]nm,wavelengths}]*)
(*][[1]];*)
(*If[NumericQ[timing],Print[If[timing>3600,ToString[timing/3600]<>"\[ThinSpace]h",IntegerString[Round[timing]]<>"\[ThinSpace]s"]]];*)


(* ::Subsection::Closed:: *)
(*maltese cross as a function of tilt angles \[Theta], \[Phi],\[CapitalDelta]n*)


(* ::Subsubsection:: *)
(*multiple runs*)


(* ::Input:: *)
(*outputDirName="dropletTilt_thetaPhi_planar_wave_ic_wvl_"<>ToString[Round[\[Lambda] 10^9]]<>"_DnLC_scan";*)
(*outputFileFormats={"h5"(*,"mat"*)};*)
(*plotQ=True;*)
(*saveDataQ=True;*)
(**)
(*(** parameter scan **)*)
(*(*{nThetas,nPhis}={200,100};*)*)
(*{nThetas,nPhis}={40,20};*)
(*tiltAnglesTheta=N@Subdivide[0,\[Pi]/2,nThetas-1];*)
(*tiltAnglesPhi=N@Subdivide[0,\[Pi]/4,nPhis-1];*)
(**)
(*(** manual specification of parameters **)*)
(*{tiltAnglesTheta,tiltAnglesPhi}=N@{{0},{0}};*)
(*{nThetas,nPhis}=Length/@{tiltAnglesTheta,tiltAnglesPhi};*)
(**)
(*\[CapitalDelta]nLCs=Subdivide[0,1.0,10-1];*)
(**)
(**)
(*(** progress indicator **)*)
(*pi=0;Print[Overlay[{ProgressIndicator[Dynamic[pi],{0,1}],Dynamic[ToString[DecimalForm[pi*100.0,{4,1}]]<>"\[ThinSpace]%"]},Alignment->"Center"]];*)
(*counter=0;*)
(*If[$VersionNumber<13,outputFileFormats=DeleteCases[outputFileFormats,"mat"]];*)
(**)
(*(** create data **)*)
(*timing=AbsoluteTiming[*)
(*intensityTiltData=Table[*)
(*setRefractionIndices[\[Lambda],\[CapitalDelta]nLC];*)
(*{xVals,yVals,zVals,nxJones,nyJones,d}=setJonesCalculusDiscretization[rDropletDim];*)
(**)
(*rotaMatrix=RotationMatrix[{{0,0,1},tiltDirector[tiltAngleTheta,tiltAnglePhi]}];*)
(*rotaMatrixInv=RotationMatrix[{tiltDirector[tiltAngleTheta,tiltAnglePhi],{0,0,1}}];*)
(**)
(*(** LC director projected on xy plane angle = marker for liquid type (LC,FC,water) **)*)
(*(** calculate extraordinary refractive index Subscript[n, e]' for LC director (projected length normal to wave vector k) **)*)
(*directorAnglesRefIndices3d=Table[*)
(*position={x,y,z};*)
(*r=Norm[position];*)
(*{directorAngle,refIndexE}=If[r<=0.99,*)
(*(** inside droplet **)*)
(*If[position . tiltDirector[tiltAngleTheta,tiltAnglePhi]>0.0,*)
(*rotatedPosition=rotaMatrixInv . position;*)
(*rotatedDirector=rotaMatrix . director[rotatedPosition];*)
(*lcAngleTheta=ArcCos[rotatedDirector[[3]]];                                           (** how much is the director tilted away from Subscript[e, z]? (\[Theta]) **)*)
(*lcAnglePhi=ArcTan[rotatedDirector[[1]],rotatedDirector[[2]]];   (** how much is the director turned in \[Phi]? **)*)
(*nLCeMod=calcXYnLCe[lcAngleTheta,\[CapitalDelta]nLC-nLCo,nLCo];*)
(*(*nLCeMod=calcXYnLCeShadow[lcAngleTheta,nLCe,nLCo];*)*)
(*{lcAnglePhi,nLCeMod},{noAngle,nOil}]*)
(*(** outside droplet **)*)
(*,{noAngle,nWater}]*)
(*,{x,xVals},{y,yVals},{z,zVals}];*)
(**)
(*(** calculate polarization state of light in LC part + FC part **)*)
(*afterDropletEField=Table[*)
(*{x,y}={xVals[[xIndex]],yVals[[yIndex]]};   (** dimless **) *)
(*totalJonesMatrix=If[x^2+y^2<1,*)
(*(** ray intersection with droplet **)*)
(*(** TODO: precompute ray intersections to only use Jones matrix in LC **)*)
(*(** 1. Jones matrix for top water **)*)
(*(*dWater=1.0-Sqrt[1.0-x^2-y^2];*)
(*singleJonesMatrixH2O=jonesMatrixH2O[dWater];*)*)
(*(** 2. Jones matrix for LC part **)*)
(*(** 3. Jones matrix for FC part; internal interface plane: n.r=zOffset **)*)
(*(*n=tiltDirector[tiltAngleTheta,tiltAnglePhi];*)
(*dFC=-(n\[LeftDoubleBracket]1\[RightDoubleBracket]x+n\[LeftDoubleBracket]2\[RightDoubleBracket]y)/n\[LeftDoubleBracket]3\[RightDoubleBracket]-Sqrt[1.0-x^2-y^2];*)
(*singleJonesMatrixFC=jonesMatrixFC[dFC];*)*)
(*(** 4. Jones matrix for bottom water, same as in 1 due to symmetry of droplet **)*)
(*directorAnglesRefIndices=Reverse@directorAnglesRefIndices3d[[xIndex,yIndex]];*)
(*singleJonesMatricesLC=Map[( *)
(*{directorAngle,refIndex}=#;*)
(*If[-\[Pi]<=directorAngle<=\[Pi],*)
(*jonesMatrixLC[directorAngle,\[Delta]LC[refIndex-nLCo],\[Delta]Iso[nLCo,d]],*)
(*jonesMatrixIsotropic[\[Delta]Iso[refIndex,d]]*)
(*]*)
(*)&,directorAnglesRefIndices];*)
(*Dot@@singleJonesMatricesLC  (** output **)*)
(*,*)
(*(** no ray intersection with droplet: only water **)*)
(*jonesMatrixH2O[2rDropletDim]   (** output **)*)
(*];*)
(*beforeDropletState=jonesMatrixHorizontalPolarizer . {1,1} initialPhaseFactor[x rDropletDim,y rDropletDim];*)
(*afterDropletState=totalJonesMatrix . beforeDropletState;*)
(*jonesMatrixVerticalPolarizer . afterDropletState*)
(*(*jonesMatrixPolarizer[0.95\[Pi]/2].afterDropletState*)*)
(*,{xIndex,nxJones}*)
(*,{yIndex,nyJones}*)
(*];*)
(*(*{phaseX,phaseY}=Arg/@afterDropletEField;*)*)
(*intensityXY=Map[Abs[#[[1]]]^2+Abs[#[[2]]]^2&,afterDropletEField,{2}]; (** effectively only intensityY due to 2nd polarizer **)*)
(*avgIntensity=Mean[Flatten[intensityXY]];*)
(**)
(**)
(*parameterString="_DnLC_"<>ToString[DecimalForm[\[CapitalDelta]nLC,{3,2}]]*)
(*<>"_wvl_"<>ToString[Round[\[Lambda] 10^9]]*)
(*<>"_theta_"<>ToString[DecimalForm[tiltAngleTheta,{4,3}]]*)
(*<>"_phi_"<>ToString[DecimalForm[tiltAnglePhi,{4,3}]];*)
(**)
(*If[plotQ,*)
(*plot=ListDensityPlot[intensityXY,InterpolationOrder->0,PlotRangePadding->None,PlotRange->All,PlotLegends->BarLegend[Automatic,LegendLabel->"intensity"],FrameLabel->{"x","y"},PlotLabel->"\[Theta] = "<>ToString[DecimalForm[tiltAngleTheta 180/\[Pi],{4,1}]]<>"\[ThinSpace]\[Degree]"*)
(*,ColorFunction->GrayLevel];*)
(*Export[FileNameJoin[{NotebookDirectory[],outputDirName,"afterDropletIntensityXY"<>parameterString<>".png"}],plot];*)
(*];*)
(**)
(*(** save data **)*)
(*If[saveDataQ,*)
(*finalEfieldState=Map[ReIm,afterDropletEField];*)
(*Export[FileNameJoin[{NotebookDirectory[],outputDirName,"afterDropletEField"<>parameterString<>"."<>#}],finalEfieldState]&/@outputFileFormats;*)
(*];*)
(**)
(*(** update progress indicator **)*)
(*pi=N[++counter/(Length[tiltAnglesTheta]Length[tiltAnglesPhi] Length[\[CapitalDelta]nLCs])];*)
(**)
(*(** output **)*)
(*{tiltAngleTheta,tiltAnglePhi,avgIntensity,afterDropletEField,intensityXY}*)
(*,{tiltAngleTheta,tiltAnglesTheta}*)
(*,{tiltAnglePhi,tiltAnglesPhi}*)
(*,{\[CapitalDelta]nLC,\[CapitalDelta]nLCs}];*)
(*][[1]];*)
(*If[NumberQ[timing],Print[If[timing>3600,ToString[timing/3600]<>"\[ThinSpace]h",ToString[timing]<>"\[ThinSpace]s"]]];*)


(* ::Subsection::Closed:: *)
(*malteser cross as a function of tilt angles \[Theta], \[Phi],\[CapitalDelta]n - raytracer**)


(* ::Subsubsection:: *)
(*multiple runs*)


(* ::Input:: *)
(*outputDirName="dropletTilt_thetaPhi_raytraced_planar_ic_wvl_"<>ToString[Round[\[Lambda] 10^9]]<>"_DnLC_scan";*)
(*outputFileFormats={"h5"(*,"mat"*)};*)
(*plotQ=True;*)
(*saveDataQ=True;*)
(**)
(*(** parameter scan **)*)
(*(*{nThetas,nPhis}={200,100};*)*)
(*{nThetas,nPhis}={40,20};*)
(*tiltAnglesTheta=N@Subdivide[0,\[Pi]/2,nThetas-1];*)
(*tiltAnglesPhi=N@Subdivide[0,\[Pi]/4,nPhis-1];*)
(**)
(*(** manual specification of parameters **)*)
(*{tiltAnglesTheta,tiltAnglesPhi}=N@{{0},{0}};*)
(*{nThetas,nPhis}=Length/@{tiltAnglesTheta,tiltAnglesPhi};*)
(**)
(*\[CapitalDelta]nLCs=Subdivide[0,1.0,1];*)
(**)
(**)
(*(** progress indicator **)*)
(*pi=0;Print[Overlay[{ProgressIndicator[Dynamic[pi],{0,1}],Dynamic[ToString[DecimalForm[pi*100.0,{4,1}]]<>"\[ThinSpace]%"]},Alignment->"Center"]];*)
(*counter=0;*)
(*If[$VersionNumber<13,outputFileFormats=DeleteCases[outputFileFormats,"mat"]];*)
(**)
(*(** create data **)*)
(*timing=AbsoluteTiming[*)
(*intensityTiltData=Table[*)
(*setRefractionIndices[\[Lambda],\[CapitalDelta]nLC];*)
(*{xVals,yVals,zVals,nxJones,nyJones,d}=setJonesCalculusDiscretization[rDropletDim];*)
(**)
(*rotaMatrix=RotationMatrix[{{0,0,1},tiltDirector[tiltAngleTheta,tiltAnglePhi]}];*)
(*rotaMatrixInv=RotationMatrix[{tiltDirector[tiltAngleTheta,tiltAnglePhi],{0,0,1}}];*)
(**)
(*(** calculate polarization state of light in LC part + FC part **)*)
(*n=tiltDirector[tiltAngleTheta,tiltAnglePhi];*)
(*afterDropletEField=Table[*)
(*{x,y}={xVals[[xIndex]],yVals[[yIndex]]};   (** dimless **) *)
(*zTop=Sqrt[1.0-x^2-y^2];*)
(*zInterface=-(n[[1]]x+n[[2]]y)/n[[3]];*)
(*zBottom=-zTop;  (** due to symmetry of droplet **)*)
(*totalJonesMatrix=If[x^2+y^2<1,*)
(*(** ray intersection with droplet **)*)
(*(** TODO: precompute ray intersections to only use Jones matrix in LC **)*)
(*(** 1. Jones matrix for top water **)*)
(*dWater=1.0-zTop;*)
(*singleJonesMatrixH2O=jonesMatrixH2O[dWater];*)
(*(** 2. Jones matrix for LC part **)*)
(*(** a) get a set of positions in LC part **)*)
(*zValsXY=Select[zVals,zInterface<#<zTop&];*)
(*dVoxel=d/rDropletDim; (** dimless **)*)
(*(** dVoxel at boundary: half distance to next regular point; voxel is centered around sample point **)*)
(*dVoxels=ConstantArray[dVoxel,Length[zValsXY]];*)
(*If[#<0.5dVoxel,*)
(*dVoxels[[1]]=#+0.5dVoxel;,*)
(*dTop=#-0.5dVoxel;PrependTo[zValsXY,zTop-0.5dTop];PrependTo[dVoxels,dTop];]&@(zTop-zValsXY[[1]]);*)
(*If[#<0.5dVoxel,*)
(*dVoxels[[-1]]=#+0.5dVoxel;,*)
(*dBottom=#-0.5dVoxel;*)
(*AppendTo[zValsXY,zInterface+0.5dBottom];*)
(*AppendTo[dVoxels,dBottom];]&@(zValsXY[[-1]]-zInterface);*)
(**)
(*positions={x,y,#}&/@zValsXY;*)
(**)
(*singleJonesMatrixLC=Dot@@MapThread[( *)
(*{directorAngle,refIndex}=getLC\[Phi]n[#1];*)
(*jonesMatrixLC[directorAngle,\[Delta]LC[refIndex-nLCo,#2],\[Delta]Iso[nLCo,#2]]*)
(*)&,{positions,dVoxels}];*)
(*(** 3. Jones matrix for FC part; internal interface plane: n.r=zOffset **)*)
(*dFC=zInterface-zBottom;*)
(*singleJonesMatrixFC=jonesMatrixFC[dFC];*)
(*(** 4. Jones matrix for bottom water, same as in 1 due to symmetry of droplet **)*)
(*Dot@@{singleJonesMatrixH2O,singleJonesMatrixFC,singleJonesMatrixLC,singleJonesMatrixH2O}  (** output **)*)
(*,*)
(*(** no ray intersection with droplet: only water **)*)
(*jonesMatrixH2O[2rDropletDim]   (** output **)*)
(*];*)
(*beforeDropletState=jonesMatrixHorizontalPolarizer . {1,1} initialPhaseFactor[x rDropletDim,y rDropletDim];*)
(*afterDropletState=totalJonesMatrix . beforeDropletState;*)
(*jonesMatrixVerticalPolarizer . afterDropletState*)
(*(*jonesMatrixPolarizer[0.95\[Pi]/2].afterDropletState*)*)
(*,{xIndex,nxJones}*)
(*,{yIndex,nyJones}*)
(*];*)
(*(*{phaseX,phaseY}=Arg/@afterDropletEField;*)*)
(*intensityXY=Map[Abs[#[[1]]]^2+Abs[#[[2]]]^2&,afterDropletEField,{2}]; (** effectively only intensityY due to 2nd polarizer **)*)
(*avgIntensity=Mean[Flatten[intensityXY]];*)
(**)
(*parameterString="_DnLC_"<>ToString[DecimalForm[\[CapitalDelta]nLC,{3,2}]]*)
(*<>"_wvl_"<>ToString[Round[\[Lambda] 10^9]]*)
(*<>"_theta_"<>ToString[DecimalForm[tiltAngleTheta,{4,3}]]*)
(*<>"_phi_"<>ToString[DecimalForm[tiltAnglePhi,{4,3}]];*)
(**)
(*If[plotQ,*)
(*plot=ListDensityPlot[intensityXY,InterpolationOrder->0,PlotRangePadding->None,PlotRange->All,PlotLegends->BarLegend[Automatic,LegendLabel->"intensity"],FrameLabel->{"x","y"},PlotLabel->"\[Theta] = "<>ToString[DecimalForm[tiltAngleTheta 180/\[Pi],{4,1}]]<>"\[ThinSpace]\[Degree]"*)
(*,ColorFunction->GrayLevel];*)
(*Export[FileNameJoin[{NotebookDirectory[],outputDirName,"afterDropletIntensityXY"<>parameterString<>".png"}],plot];*)
(*];*)
(**)
(*(** save data **)*)
(*If[saveDataQ,*)
(*finalEfieldState=Map[ReIm,afterDropletEField];*)
(*Export[FileNameJoin[{NotebookDirectory[],outputDirName,"afterDropletEField"<>parameterString<>"."<>#}],finalEfieldState]&/@outputFileFormats;*)
(*];*)
(**)
(*(** update progress indicator **)*)
(*pi=N[++counter/(Length[tiltAnglesTheta]Length[tiltAnglesPhi] Length[\[CapitalDelta]nLCs])];*)
(**)
(*(** output **)*)
(*{tiltAngleTheta,tiltAnglePhi,avgIntensity,afterDropletEField,intensityXY}*)
(*,{tiltAngleTheta,tiltAnglesTheta}*)
(*,{tiltAnglePhi,tiltAnglesPhi}*)
(*,{\[CapitalDelta]nLC,\[CapitalDelta]nLCs}];*)
(*][[1]];*)
(*If[NumberQ[timing],Print[If[timing>3600,ToString[timing/3600]<>"\[ThinSpace]h",ToString[timing]<>"\[ThinSpace]s"]]];*)


(* ::Subsection::Closed:: *)
(*malteser cross as a function of tilt angles \[Theta], \[Phi], rotating cross polarizers*)


(* ::Subsubsection:: *)
(*main*)


(* ::Input:: *)
(*(*tiltAnglesTheta={0.3};*)
(*polarizerAnglesPhi=N@{\[Pi]/8};*)*)
(*tiltAnglePhi=0.0;*)
(*(*tiltAnglesTheta=N@Subdivide[0,\[Pi]/2,40];*)
(*polarizerAnglesPhi=N@Subdivide[0,\[Pi]/4,20];*)*)
(*tiltAnglesTheta=N@{\[Pi]/4};*)
(*polarizerAnglesPhi=N@{\[Pi]/4};*)
(**)
(*Print[AbsoluteTiming[*)
(*intensityTiltDataPolarizer=Table[*)
(*rotaMatrix=RotationMatrix[{{0,0,1},tiltDirector[tiltAngleTheta,tiltAnglePhi]}];*)
(*rotaMatrixInv=RotationMatrix[{tiltDirector[tiltAngleTheta,tiltAnglePhi],{0,0,1}}];*)
(*rotaMatrixPhi=RotationMatrix[polarizerAnglePhi];*)
(*rotaMatrixPhiInv=RotationMatrix[-polarizerAnglePhi];*)
(**)
(*(** LC director projected on xy plane angle = marker for liquid type (LC,FC,water) **)*)
(*(** calculate extraordinary refractive index Subscript[n, e]' for LC director (projected length normal to wave vector k) **)*)
(*directorAnglesRefIndices3d=Table[*)
(*position={x,y,z};*)
(*{directorAngle,refIndex}=If[Norm[position]<=0.99,*)
(*If[position . tiltDirector[tiltAngleTheta,tiltAnglePhi]>0.0,*)
(*rotatedPosition=rotaMatrixInv . position;*)
(*rotatedDirector=rotaMatrix . director[rotatedPosition];*)
(*lcAngle=ArcTan[rotatedDirector[[1]],rotatedDirector[[2]]];*)
(*nLCeMod=calcXYProjectedDirectorLength[rotatedDirector,waveVector,nLCe];*)
(*{lcAngle,nLCeMod},{noAngle,nOil}],{noAngle,nWater}]*)
(*,{x,xVals},{y,yVals},{z,zVals}];*)
(**)
(*Table[*)
(*(** calculate polarization state of light in LC part + FC part **)*)
(*afterDropletEField=Table[*)
(*directorAnglesRefIndices=Reverse@directorAnglesRefIndices3d[[xIndex,yIndex]];*)
(*(** ray intersection with the droplet **)*)
(*singleJonesMatrices=Map[( *)
(*{directorAngle,refIndex}=#;*)
(*If[-\[Pi]<=directorAngle<=\[Pi],*)
(*jonesMatrixLC[directorAngle,\[Delta]LC[refIndex-nLCo],\[Delta]Iso[nLCo]],*)
(*jonesMatrixIsotropic[\[Delta]Uni[refIndex]]*)
(*]*)
(*)&,directorAnglesRefIndices];*)
(*dropletJonesMatrix=Dot@@singleJonesMatrices;*)
(*afterDropletState=dropletJonesMatrix . (rotaMatrixPhi . polarizationStateHorizontal);*)
(*(rotaMatrixPhi . jonesMatrixVerticalPolarizer . rotaMatrixPhiInv) . afterDropletState*)
(*,{xIndex,nx}*)
(*,{yIndex,ny}*)
(*];*)
(*intensity=Map[Abs[#[[1]]]^2+Abs[#[[2]]]^2&,afterDropletEField,{2}];*)
(*avgIntensity=Mean[Flatten[intensity]];*)
(**)
(*(** save data **)*)
(*finalEfieldState=Map[ReIm,afterDropletEField];*)
(*Export[FileNameJoin[{NotebookDirectory[],"dropletTilt_thetaPhi_lcfch2o_4_rotatingPolarizer","afterDropletEField_550_theta_"<>ToString[DecimalForm[tiltAngleTheta,{4,3}]]<>"_phi_"<>ToString[DecimalForm[polarizerAnglePhi,{4,3}]]<>".h5"}],finalEfieldState];*)
(**)
(*(** output **)*)
(*{tiltAngleTheta,polarizerAnglePhi,avgIntensity,afterDropletEField,intensity}*)
(*,{polarizerAnglePhi,polarizerAnglesPhi}]*)
(*,{tiltAngleTheta,tiltAnglesTheta}];*)
(*][[1]],"\[ThinSpace]s"];*)


(* ::Text:: *)
(*previously : *)
(*11575.8042653`"\[ThinSpace]s"*)
(*38239.0997881` s*)


(* ::Subsubsection::Closed:: *)
(*save E fields*)


(* ::Input:: *)
(*(** save output electric field **)*)
(*saveStateQ=True;*)
(*If[saveStateQ,*)
(*Do[*)
(*finalEfieldState=Map[ReIm,intensityTiltData[[i,j,4]]];*)
(*Export[FileNameJoin[{$HomeDirectory,"Desktop","dropletTiltTest_brightSpot","afterDropletEField_550_i_"<>ToString[i]<>"_j_"<>ToString[j]<>".mat"}],finalEfieldState];*)
(*,{i,Length[tiltAnglesTheta]},{j,Length[tiltAnglesTheta]}]*)
(*];*)


(* ::Subsubsection::Closed:: *)
(*show intensity*)


(* ::Input:: *)
(*i=-1;*)
(*plot=ListDensityPlot[intensityTiltData[[1,1,5]],InterpolationOrder->0,PlotRangePadding->None,PlotRange->All,PlotLegends->BarLegend[Automatic,LegendLabel->"intensity"],FrameLabel->{"x","y"},PlotLabel->"\[Theta] = "<>ToString[DecimalForm[intensityTiltData[[i,1,1]]180/\[Pi],{4,1}]]<>" \[Degree]"*)
(*,ColorFunction->GrayLevel*)
(*(*,PerformanceGoal->"Speed"*)]*)


(* ::Subsubsection::Closed:: *)
(*scalar intensities right after droplet - lookup map*)


(* ::Input:: *)
(*plot=Show[*)
(*ListPlot3D[{180/\[Pi]#[[1]],180/\[Pi]#[[2]],#[[3]]}&/@Flatten[intensityTiltData[[All,All,{1,2,3}]],1],AxesLabel->{"\[Theta] (\[Degree])","\[Phi] (\[Degree])","I"},AxesStyle->Black,Mesh->None],*)
(*ListPointPlot3D[{180/\[Pi]#[[1]],180/\[Pi]#[[2]],#[[3]]}&/@Flatten[intensityTiltData[[All,All,{1,2,3}]],1],AxesLabel->{"\[Theta] (\[Degree])","\[Phi] (\[Degree])","I"},AxesStyle->Black,PlotStyle->Black]*)
(*]*)
(*Export[FileNameJoin[{NotebookDirectory[],"thetaPhi_intensity.png"}],plot]*)


(* ::Subsubsection::Closed:: *)
(*create movie of tilted intensities for different resolutions*)


(* ::Input:: *)
(*Manipulate[*)
(*ListDensityPlot[intensityTiltData[[i,j,5]],InterpolationOrder->0,PlotRangePadding->None,PlotRange->All,PlotLegends->BarLegend[Automatic,LegendLabel->"intensity"],FrameLabel->{"x","y"},PlotLabel->"\[Theta] = "<>ToString@intensityTiltData[[i,j,1]],ColorFunction->GrayLevel]*)
(*,{i,1,Length[tiltAnglesTheta],1}*)
(*,{j,1,Length[tiltAnglesPhi],1}*)
(*]*)


(* ::Input:: *)
(*counter=0;*)
(*Do[*)
(*label="\[Theta] = "<>ToString[DecimalForm[intensityTiltData[[i,1]]180/\[Pi],{4,1}]]<>"\[Degree], \[Phi] = "<>ToString[DecimalForm[tiltAnglePhi 180/\[Pi],{4,1}]]<>"\[Degree]";*)
(*plotHighRes=ListDensityPlot[intensityTiltData[[i,j,4]],InterpolationOrder->0,PlotRangePadding->None,PlotRange->All,PlotLegends->BarLegend[Automatic,LegendLabel->"intensity"],FrameLabel->{"x","y"},ColorFunction->GrayLevel,ImageSize->{Automatic,400}];*)
(*plotLowRes=ListDensityPlot[intensityTiltData[[i,j,4]],InterpolationOrder->0,PlotRangePadding->None,PlotRange->All,PlotLegends->BarLegend[Automatic,LegendLabel->"intensity"],FrameLabel->{"x","y"},ColorFunction->GrayLevel,ImageSize->{Automatic,400},PerformanceGoal->"Speed"];*)
(*plot=Grid[{{label,SpanFromLeft},{plotHighRes,plotLowRes}}];*)
(*Export[FileNameJoin[{NotebookDirectory[],"mma_polarizationMovie_"<>ToString[tiltAnglePhi 180/\[Pi]]<>"_deg","p_"<>IntegerString[++counter,10,IntegerLength[Length@intensityTiltData]+1]<>".png"}],plot];*)
(*,{i,Length[intensityTiltData]}];*)


(* ::Subsubsection::Closed:: *)
(*width of composite Lorentzian*)


(* ::Input:: *)
(*Plot[{1/(1+x^2),1/(1+(x/2.0)^2),1/(1+(x/1.5)^2),0.5(1/(1+x^2)+1/(1+(x/2.0)^2))},{x,0,10},PlotRange->All]*)


(* ::Subsubsection::Closed:: *)
(*debug*)


(* ::Input:: *)
(*avgRefIndexE=Table[*)
(*Mean[directorAnglesRefIndices3d[[xIndex,yIndex,All,2]]]*)
(*,{xIndex,nx},{yIndex,ny}*)
(*];*)


(* ::Input:: *)
(*ListDensityPlot[avgRefIndexE,InterpolationOrder->0,PlotRange->All,PlotLegends->Automatic]*)


(* ::Input:: *)
(*ListPlot[avgRefIndexE[[50]]]*)


(* ::Input:: *)
(*totalRefIndexE=Table[*)
(*Mod[Total[directorAnglesRefIndices3d[[xIndex,yIndex,All,2]]]dz rDropletDim 10^6,\[Lambda] 10^6]*)
(*,{xIndex,nx},{yIndex,ny}*)
(*];*)
(*ListDensityPlot[totalRefIndexE,InterpolationOrder->0,PlotRange->All,PlotLegends->Automatic]*)


(* ::Input:: *)
(*Histogram[Flatten[directorAnglesRefIndices3d[[All,All,All,3]]]]*)


(* ::Input:: *)
(*\[CapitalDelta]n=directorAnglesRefIndices3d[[50,All,All,2]]\[Transpose]-nLCo;*)


(* ::Input:: *)
(*ListDensityPlot[\[CapitalDelta]n,PlotLegends->Automatic,PlotRange->All,InterpolationOrder->0]*)


(* ::Input:: *)
(*avgnLC=Table[*)
(*Mean[directorAnglesRefIndices3d[[xIndex,yIndex,All,2]]]*)
(*,{xIndex,Dimensions[directorAnglesRefIndices3d][[1]]}*)
(*,{yIndex,Dimensions[directorAnglesRefIndices3d][[2]]}];*)
(*ListDensityPlot[avgnLC,InterpolationOrder->0]*)


(* ::Input:: *)
(*ListContourPlot[avgnLC,Contours->{nz nWater,1400,1410(*nz(0.5 nLCo+0.5nLCe)*)}]*)


(* ::Subsection::Closed:: *)
(*malteser cross as a function of tilt angles \[Theta], \[Phi] (via optical path length)**)


(* ::Text:: *)
(*old runtime: *)


(* ::Subsubsection:: *)
(*main*)


(* ::Input:: *)
(*(*{nThetas,nPhis}={200,100};*)
(*{nThetas,nPhis}={40,20};*)*)
(*{nThetas,nPhis}={2,2};*)
(**)
(*(*tiltAnglesTheta=N[Subdivide[0,\[Pi]/2,20]];*)
(*tiltAnglesPhi=N[Subdivide[0,\[Pi]/4,10]];*)*)
(*nThetas=nPhis=1;*)
(*(*tiltAnglesTheta={0.0};*)
(*tiltAnglesPhi={0.0};*)*)
(*tiltAnglesTheta=N@{\[Pi]/4};*)
(*tiltAnglesPhi=N@{\[Pi]/4};*)
(*(*tiltAnglesTheta=N@{\[Pi]/32,\[Pi]/16,\[Pi]/8};*)
(*tiltAnglesPhi=N@{0,\[Pi]/8,\[Pi]/4};*)*)
(*(*tiltAnglesTheta=N@Subdivide[0,\[Pi]/2,nThetas-1];*)
(*tiltAnglesPhi=N@Subdivide[0,\[Pi]/4,nPhis-1];*)*)
(**)
(*(** progress indicator **)*)
(*pi=0;Print[Overlay[{ProgressIndicator[Dynamic[pi],{0,1}],Dynamic[ToString[DecimalForm[pi*100,{4,1}]]<>"\[ThinSpace]%"]},Alignment->"Center"]];*)
(*counter=0;*)
(**)
(*(** create data **)*)
(*Print[AbsoluteTiming[*)
(*intensityTiltData=Table[*)
(*rotaMatrix=RotationMatrix[{{0,0,1},tiltDirector[tiltAngleTheta,tiltAnglePhi]}];*)
(*rotaMatrixInv=RotationMatrix[{tiltDirector[tiltAngleTheta,tiltAnglePhi],{0,0,1}}];*)
(**)
(*(** LC director projected on xy plane angle = marker for liquid type (LC,FC,water) **)*)
(*(** calculate extraordinary refractive index Subscript[n, e]' for LC director (projected length normal to wave vector k) **)*)
(*directorAnglesRefIndices3d=Table[*)
(*position={x,y,z};*)
(*r=Norm[position];*)
(*{directorAngle,refIndexE}=If[r<=0.99,*)
(*(** inside droplet **)*)
(*If[position . tiltDirector[tiltAngleTheta,tiltAnglePhi]>0.0,*)
(*rotatedPosition=rotaMatrixInv . position;*)
(*rotatedDirector=rotaMatrix . director[rotatedPosition];*)
(*lcAngle=ArcTan[rotatedDirector[[1]],rotatedDirector[[2]]];*)
(*nLCeMod=calcXYProjectedDirectorLength[rotatedDirector,waveVector,nLCe];*)
(*{lcAngle,nLCeMod},{noAngle,nOil}]*)
(*(** outside droplet **)*)
(*,{noAngle,nWater}]*)
(*,{x,xVals},{y,yVals},{z,zVals}];*)
(**)
(*(** calculate polarization state of light in LC part + FC part **)*)
(*afterDropletEField=Table[*)
(*directorAnglesRefIndices=Reverse@directorAnglesRefIndices3d[[xIndex,yIndex]];*)
(*(** ray intersection with the droplet **)*)
(*singleJonesMatrices=Map[( *)
(*{directorAngle,refIndex}=#;*)
(*If[-\[Pi]<=directorAngle<=\[Pi],*)
(*jonesMatrixLC[directorAngle,\[Delta]LC[refIndex-nLCo],\[Delta]Iso[nLCo]],*)
(*jonesMatrixIsotropic[\[Delta]Uni[refIndex]]*)
(*]*)
(*)&,directorAnglesRefIndices];*)
(*dropletJonesMatrix=Dot@@singleJonesMatrices;*)
(*afterDropletState=dropletJonesMatrix . polarizationStateHorizontal;*)
(*jonesMatrixVerticalPolarizer . afterDropletState*)
(*,{xIndex,nxJones}*)
(*,{yIndex,nyJones}*)
(*];*)
(*(*{phaseX,phaseY}=Arg/@afterDropletEField;*)*)
(*intensity=Map[Abs[#[[1]]]^2+Abs[#[[2]]]^2&,afterDropletEField,{2}];*)
(*avgIntensity=Mean[Flatten[intensity]];*)
(**)
(*(** save data **)*)
(*finalEfieldState=Map[ReIm,afterDropletEField];*)
(*pthout=FileNameJoin[{NotebookDirectory[],"dropletTilt_thetaPhi_lcfch2o_4"}];*)
(*Export[FileNameJoin[{pthout,"afterDropletEField_550_theta_"<>ToString[DecimalForm[tiltAngleTheta,{4,3}]]<>"_phi_"<>ToString[DecimalForm[tiltAnglePhi,{4,3}]]<>".h5"}],finalEfieldState];*)
(**)
(*(** update progress indicator **)*)
(*pi=N[++counter/(Length[tiltAnglesTheta]Length[tiltAnglesPhi])];*)
(**)
(*(** output **)*)
(*{tiltAngleTheta,tiltAnglePhi,avgIntensity,afterDropletEField,intensity}*)
(*,{tiltAngleTheta,tiltAnglesTheta}*)
(*,{tiltAnglePhi,tiltAnglesPhi}];*)
(*][[1]]];*)


(* ::Input:: *)
(*38239.0997881`*)


(* ::Subsubsection::Closed:: *)
(*save E fields*)


(* ::Input:: *)
(*(** save output electric field **)*)
(*saveStateQ=True;*)
(*If[saveStateQ,*)
(*Do[*)
(*finalEfieldState=Map[ReIm,intensityTiltData[[i,j,4]]];*)
(*Export[FileNameJoin[{$HomeDirectory,"Desktop","dropletTiltTest_brightSpot","afterDropletEField_550_i_"<>ToString[i]<>"_j_"<>ToString[j]<>".h5"}],finalEfieldState];*)
(*,{i,Length[tiltAnglesTheta]},{j,Length[tiltAnglesTheta]}]*)
(*];*)


(* ::Subsubsection:: *)
(*show intensity*)


(* ::Input:: *)
(*i=-1;*)
(*plot=ListDensityPlot[intensityTiltDataPolarizer[[1,1,5]],InterpolationOrder->0,PlotRangePadding->None,PlotRange->All,PlotLegends->BarLegend[Automatic,LegendLabel->"intensity"],FrameLabel->{"x","y"},PlotLabel->"\[Theta] = "<>ToString[DecimalForm[intensityTiltData[[i,1,1]]180/\[Pi],{4,1}]]<>" \[Degree]"*)
(*,ColorFunction->GrayLevel*)
(*(*,PerformanceGoal->"Speed"*)]*)


(* ::Subsubsection::Closed:: *)
(*scalar intensities right after droplet - lookup map*)


(* ::Input:: *)
(*plot=Show[*)
(*ListPlot3D[{180/\[Pi]#[[1]],180/\[Pi]#[[2]],#[[3]]}&/@Flatten[intensityTiltData[[All,All,{1,2,3}]],1],AxesLabel->{"\[Theta] (\[Degree])","\[Phi] (\[Degree])","I"},AxesStyle->Black,Mesh->None],*)
(*ListPointPlot3D[{180/\[Pi]#[[1]],180/\[Pi]#[[2]],#[[3]]}&/@Flatten[intensityTiltData[[All,All,{1,2,3}]],1],AxesLabel->{"\[Theta] (\[Degree])","\[Phi] (\[Degree])","I"},AxesStyle->Black,PlotStyle->Black]*)
(*]*)
(*Export[FileNameJoin[{NotebookDirectory[],"thetaPhi_intensity.png"}],plot]*)


(* ::Subsubsection::Closed:: *)
(*create movie of tilted intensities for different resolutions*)


(* ::Input:: *)
(*Manipulate[*)
(*ListDensityPlot[intensityTiltData[[i,j,5]],InterpolationOrder->0,PlotRangePadding->None,PlotRange->All,PlotLegends->BarLegend[Automatic,LegendLabel->"intensity"],FrameLabel->{"x","y"},PlotLabel->"\[Theta] = "<>ToString@intensityTiltData[[i,j,1]],ColorFunction->GrayLevel]*)
(*,{i,1,Length[tiltAnglesTheta],1}*)
(*,{j,1,Length[tiltAnglesPhi],1}*)
(*]*)


(* ::Input:: *)
(*counter=0;*)
(*Do[*)
(*label="\[Theta] = "<>ToString[DecimalForm[intensityTiltData[[i,1]]180/\[Pi],{4,1}]]<>"\[Degree], \[Phi] = "<>ToString[DecimalForm[tiltAnglePhi 180/\[Pi],{4,1}]]<>"\[Degree]";*)
(*plotHighRes=ListDensityPlot[intensityTiltData[[i,j,4]],InterpolationOrder->0,PlotRangePadding->None,PlotRange->All,PlotLegends->BarLegend[Automatic,LegendLabel->"intensity"],FrameLabel->{"x","y"},ColorFunction->GrayLevel,ImageSize->{Automatic,400}];*)
(*plotLowRes=ListDensityPlot[intensityTiltData[[i,j,4]],InterpolationOrder->0,PlotRangePadding->None,PlotRange->All,PlotLegends->BarLegend[Automatic,LegendLabel->"intensity"],FrameLabel->{"x","y"},ColorFunction->GrayLevel,ImageSize->{Automatic,400},PerformanceGoal->"Speed"];*)
(*plot=Grid[{{label,SpanFromLeft},{plotHighRes,plotLowRes}}];*)
(*Export[FileNameJoin[{NotebookDirectory[],"mma_polarizationMovie_"<>ToString[tiltAnglePhi 180/\[Pi]]<>"_deg","p_"<>IntegerString[++counter,10,IntegerLength[Length@intensityTiltData]+1]<>".png"}],plot];*)
(*,{i,Length[intensityTiltData]}];*)


(* ::Subsubsection::Closed:: *)
(*debug*)


(* ::Input:: *)
(*avgRefIndexE=Table[*)
(*Mean[directorAnglesRefIndices3d[[xIndex,yIndex,All,2]]]*)
(*,{xIndex,nx},{yIndex,ny}*)
(*];*)


(* ::Input:: *)
(*ListDensityPlot[avgRefIndexE,InterpolationOrder->0,PlotRange->All,PlotLegends->Automatic]*)


(* ::Input:: *)
(*ListPlot[avgRefIndexE[[50]]]*)


(* ::Input:: *)
(*totalRefIndexE=Table[*)
(*Mod[Total[directorAnglesRefIndices3d[[xIndex,yIndex,All,2]]]dz rDropletDim 10^6,\[Lambda] 10^6]*)
(*,{xIndex,nx},{yIndex,ny}*)
(*];*)
(*ListDensityPlot[totalRefIndexE,InterpolationOrder->0,PlotRange->All,PlotLegends->Automatic]*)


(* ::Input:: *)
(*Histogram[Flatten[directorAnglesRefIndices3d[[All,All,All,3]]]]*)


(* ::Input:: *)
(*\[CapitalDelta]n=directorAnglesRefIndices3d[[50,All,All,2]]\[Transpose]-nLCo;*)


(* ::Input:: *)
(*ListDensityPlot[\[CapitalDelta]n,PlotLegends->Automatic,PlotRange->All,InterpolationOrder->0]*)


(* ::Input:: *)
(*avgnLC=Table[*)
(*Mean[directorAnglesRefIndices3d[[xIndex,yIndex,All,2]]]*)
(*,{xIndex,Dimensions[directorAnglesRefIndices3d][[1]]}*)
(*,{yIndex,Dimensions[directorAnglesRefIndices3d][[2]]}];*)
(*ListDensityPlot[avgnLC,InterpolationOrder->0]*)


(* ::Input:: *)
(*ListContourPlot[avgnLC,Contours->{nz nWater,1400,1410(*nz(0.5 nLCo+0.5nLCe)*)}]*)


(* ::Subsubsection::Closed:: *)
(*old (without water)*)


(* ::Input:: *)
(*(*afterDropletState=dropletJonesMatrix.RotationMatrix[polarizerAnglePhi].polarizationStateHorizontal;*)
(*Inverse[RotationMatrix[polarizerAnglePhi]].jonesMatrixVerticalPolarizer.RotationMatrix[polarizerAnglePhi].afterDropletState*)
(**)
(*directorAnglesRefIndices=Reverse@directorAnglesRefIndices3d[[xIndex,yIndex]];*)
(*filteredAnglesRefIndices=If[#\[LeftDoubleBracket]1\[RightDoubleBracket]!=waterAngle,#,Nothing]&/@directorAnglesRefIndices;*)
(*If[filteredAnglesRefIndices!={},*)
(*(** ray intersection with the droplet **)*)
(*(*singleJonesMatrices=Map[If[#\[LeftDoubleBracket]1\[RightDoubleBracket]!=fcAngle,jonesMatrixLC[#\[LeftDoubleBracket]1\[RightDoubleBracket],nLCe],jonesMatrixFC]&,filteredAnglesRefIndices];*)*)
(*singleJonesMatrices=Map[If[#\[LeftDoubleBracket]1\[RightDoubleBracket]!=fcAngle,jonesMatrixLC[#\[LeftDoubleBracket]1\[RightDoubleBracket],(*nLCe*)#\[LeftDoubleBracket]2\[RightDoubleBracket]],jonesMatrixFC]&,filteredAnglesRefIndices];dropletJonesMatrix=Dot@@singleJonesMatrices;*)
(*afterDropletState=dropletJonesMatrix.polarizationStateHorizontal;*)
(*jonesMatrixVerticalPolarizer.afterDropletState,*)
(*(** no ray intersection with the droplet **)*)
(*jonesMatrixVerticalPolarizer.polarizationStateHorizontal*)
(*]*)*)


(* ::Subsection::Closed:: *)
(*malteser cross as a function of tilt angles \[Theta], \[Phi], non-vertical rays*)


(* ::Subsubsection:: *)
(*main*)


(* ::Input:: *)
(*tiltAnglesTheta=N[Subdivide[0,\[Pi]/2,20]];*)
(*tiltAnglesPhi=N[Subdivide[0,\[Pi]/4,10]];*)
(*tiltAnglesTheta={0.0};*)
(*tiltAnglesPhi={0.0};*)
(*tiltDirector[\[Theta]_,\[Phi]_]:={Cos[\[Phi]]Sin[\[Theta]],Sin[\[Phi]]Sin[\[Theta]],Cos[\[Theta]]};*)
(*saveStateQ=False;*)
(*Print[AbsoluteTiming[*)
(*intensityTiltData=Table[*)
(*rotaMatrix=RotationMatrix[{{0,0,1},tiltDirector[tiltAngleTheta,tiltAnglePhi]}];*)
(*rotaMatrixInv=RotationMatrix[{tiltDirector[tiltAngleTheta,tiltAnglePhi],{0,0,1}}];*)
(**)
(*(** LC director projected on xy plane angle = marker for liquid type (LC,FC,water) **)*)
(*(** calculate extraordinary refractive index Subscript[n, e]'' for LC director != wave director, weird assumption: polarization does not change as light goes through droplet even though it is a wave retarder **)*)
(*directorAnglesRefIndices3d=Table[*)
(*position={x,y,z};*)
(*{directorAngle,refIndex}=If[Norm[position]<=0.99,*)
(*If[position . tiltDirector[tiltAngleTheta,tiltAnglePhi]>0.0,*)
(*rotatedPosition=rotaMatrixInv . position;*)
(*rotatedDirector=rotaMatrix . director[rotatedPosition];*)
(*lcAngle=ArcTan[rotatedDirector[[1]],rotatedDirector[[2]]];*)
(*nLC=calcLCrefIndex[rotatedDirector,polarizationVector];*)
(*{lcAngle,nLC},{fcAngle,nOil}],{waterAngle,nWater}];*)
(*{directorAngle,refIndex}*)
(*,{x,xVals},{y,yVals},{z,zVals}];*)
(**)
(*(** calculate polarization state of light in LC part + FC part **)*)
(*dropletPolarizationState=Table[*)
(*directorAnglesRefIndices=Reverse@directorAnglesRefIndices3d[[xIndex,yIndex,All,{1,2}]];*)
(*filteredAnglesRefIndices=Select[directorAnglesRefIndices,#[[1]]!=waterAngle&];*)
(*afterVerticalPolarizerState=If[filteredAnglesRefIndices!={},*)
(*(** ray intersection with the droplet **)*)
(*singleJonesMatrices=Map[If[#[[1]]!=fcAngle,jonesMatrixLC[#[[1]],(*nLCe*)#[[2]]],jonesMatrixFC]&,filteredAnglesRefIndices];*)
(*dropletJonesMatrix=Dot@@singleJonesMatrices;*)
(*afterDropletState=dropletJonesMatrix . polarizationStateHorizontal;*)
(*jonesMatrixVerticalPolarizer . afterDropletState,*)
(*(** no ray intersection with the droplet **)*)
(*jonesMatrixVerticalPolarizer . polarizationStateHorizontal*)
(*];*)
(*afterVerticalPolarizerState*)
(*,{xIndex,Dimensions[directorAnglesRefIndices3d][[1]]},{yIndex,Dimensions[directorAnglesRefIndices3d][[2]]}];*)
(*intensity=Map[Abs[#[[1]]]^2+Abs[#[[2]]]^2&,dropletPolarizationState,{2}];*)
(**)
(*(** output **)*)
(*{tiltAngleTheta,tiltAnglePhi,Mean[Flatten[intensity]],dropletPolarizationState,intensity}*)
(*,{tiltAngleTheta,tiltAnglesTheta}*)
(*,{tiltAnglePhi,tiltAnglesPhi}];*)
(*][[1]]];*)
(**)
(*(** save polarization state **)*)
(*If[saveStateQ,*)
(*Do[*)
(*finalPolarizationState=Map[ReIm,intensityTiltData[[i,j,4]]];*)
(*Export[FileNameJoin[{$HomeDirectory,"Desktop","afterDropletPolarizationState_i_"<>ToString[i]<>"_j_"<>ToString[j]<>".mat"}],finalPolarizationState];*)
(*,{i,Length[tiltAnglesTheta]},{j,Length[tiltAnglesTheta]}]*)
(*];*)


(* ::Input:: *)
(*avgnLC=Table[*)
(*Mean[directorAnglesRefIndices3d[[xIndex,yIndex,All,2]]]*)
(*,{xIndex,Dimensions[directorAnglesRefIndices3d][[1]]}*)
(*,{yIndex,Dimensions[directorAnglesRefIndices3d][[2]]}];*)
(*ListDensityPlot[avgnLC,InterpolationOrder->0]*)


(* ::Input:: *)
(*ListContourPlot[avgnLC,Contours->{nz nWater,1400,1410(*nz(0.5 nLCo+0.5nLCe)*)}]*)


(* ::Input:: *)
(*plot=Show[*)
(*Plot[0.35Cos[\[Theta]]^2,{\[Theta],0,\[Pi]/2},PlotStyle->LightGray,PlotLegends->{"\!\(\*SuperscriptBox[\(cos\), \(2\)]\)\[Theta]"}],*)
(*ListLinePlot[intensityTiltData[[All,{1,2}]],Mesh->All,PlotLegends->{"Jones matrix"}]*)
(*,Frame->True,FrameLabel->{"tilt angle \[Theta]","intensity"},PlotLabel->Style["\[Phi] = "<>ToString[Round[tiltAnglePhi 180/\[Pi]]]<>"\[Degree]",Black,20],FrameStyle->Directive[Black,18,AbsoluteThickness[2]],AspectRatio->1*)
(*]*)
(*Export[FileNameJoin[{NotebookDirectory[],"intensity_over_tilt_phi_"<>ToString[IntegerString[tiltAnglePhi 180/\[Pi]]]<>".png"}],plot]*)


(* ::Input:: *)
(*i=-1;*)
(*ListDensityPlot[intensityTiltData[[1,1,5]],InterpolationOrder->0,PlotRangePadding->None,PlotRange->All,PlotLegends->BarLegend[Automatic,LegendLabel->"intensity"],FrameLabel->{"x","y"},PlotLabel->"\[Theta] = "<>ToString[DecimalForm[intensityTiltData[[i,1,1]]180/\[Pi],{4,1}]]<>" \[Degree]"*)
(*,ColorFunction->GrayLevel*)
(*(*,PerformanceGoal->"Speed"*)]*)


(* ::Input:: *)
(*counter=0;*)
(*Do[*)
(*label="\[Theta] = "<>ToString[DecimalForm[intensityTiltData[[i,1]]180/\[Pi],{4,1}]]<>"\[Degree], \[Phi] = "<>ToString[DecimalForm[tiltAnglePhi 180/\[Pi],{4,1}]]<>"\[Degree]";*)
(*plotHighRes=ListDensityPlot[intensityTiltData[[i,3]],InterpolationOrder->0,PlotRangePadding->None,PlotRange->All,PlotLegends->BarLegend[Automatic,LegendLabel->"intensity"],FrameLabel->{"x","y"},ColorFunction->GrayLevel,ImageSize->{Automatic,400}];*)
(*plotLowRes=ListDensityPlot[intensityTiltData[[i,3]],InterpolationOrder->0,PlotRangePadding->None,PlotRange->All,PlotLegends->BarLegend[Automatic,LegendLabel->"intensity"],FrameLabel->{"x","y"},ColorFunction->GrayLevel,ImageSize->{Automatic,400},PerformanceGoal->"Speed"];*)
(*plot=Grid[{{label,SpanFromLeft},{plotHighRes,plotLowRes}}];*)
(*Export[FileNameJoin[{NotebookDirectory[],"mma_polarizationMovie_"<>ToString[tiltAnglePhi 180/\[Pi]]<>"_deg","p_"<>IntegerString[++counter,10,IntegerLength[Length@intensityTiltData]+1]<>".png"}],plot];*)
(*,{i,Length[intensityTiltData]}];*)


(* ::Input:: *)
(*Manipulate[*)
(*ListDensityPlot[intensityTiltData[[i,3]],InterpolationOrder->0,PlotRangePadding->None,PlotRange->All,PlotLegends->BarLegend[Automatic,LegendLabel->"intensity"],FrameLabel->{"x","y"},PlotLabel->"\[Theta] = "<>ToString@intensityTiltData[[i,1]],ColorFunction->GrayLevel]*)
(*,{i,1,21,1}]*)


(* ::Subsection::Closed:: *)
(*compare optical models*)


(* ::Input:: *)
(*plot=Show[*)
(*Plot[0.35Cos[\[Theta]]^2,{\[Theta],0,\[Pi]/2},PlotStyle->LightGray,PlotLegends->{"\!\(\*SuperscriptBox[\(cos\), \(2\)]\)\[Theta]"}],*)
(*ListLinePlot[intensityTiltData[[All,{1,2}]],Mesh->All,PlotLegends->{"Jones matrix"}]*)
(*,Frame->True,FrameLabel->{"tilt angle \[Theta]","intensity"},PlotLabel->Style["\[Phi] = "<>ToString[Round[tiltAnglePhi 180/\[Pi]]]<>"\[Degree]",Black,20],FrameStyle->Directive[Black,18,AbsoluteThickness[2]],AspectRatio->1*)
(*]*)
(*Export[FileNameJoin[{NotebookDirectory[],"intensity_over_tilt_phi_"<>ToString[IntegerString[tiltAnglePhi 180/\[Pi]]]<>".png"}],plot]*)


(* ::Subsection::Closed:: *)
(*calculate refractive index / optical length 3d*)


(* ::Input:: *)
(*Clear[x,y,z,director,nLC];*)
(*(*director[x_,y_,z_]=Normalize[{u[x,y,z],v[x,y,z],w[x,y,z]}/.solution];*)*)
(*(*nLCBad[director_]:=Norm[surfaceProjection[-{0,0,1},director]];*)*)
(*(*nLCBad2[cos\[Theta]_]:=1.0/Sqrt[cos\[Theta]^2/nLCo^2+(1-cos\[Theta]^2)/nLCe^2];*)*)


(* ::Input:: *)
(*{dx,dy,dz}=0.02{1,1,1};*)
(*tiltAngle=0 \[Pi](*\[Pi]/4*);*)
(**)
(*Print[AbsoluteTiming[*)
(*refractiveIndexData3d=Table[*)
(*refIndex=If[Norm[{x,y,z}]<=0.99,If[z Cos[tiltAngle]+y Sin[tiltAngle]>0,*)
(*rotatedDirector=RotationMatrix[-tiltAngle,{1,0,0}] . director[x,y Cos[tiltAngle]-z Sin[tiltAngle],z Cos[tiltAngle]+y Sin[tiltAngle]];*)
(*(*cosDirectorTilt={0,0,-1}.rotatedDirector;*)*)
(*(*nLC[cosDirectorTilt]*)*)
(*nLC[rotatedDirector,polarizationVector],nOil],nWater];*)
(*{x,y,z,refIndex}*)
(*,{x,-1,1,dx},{y,-1,1,dy},{z,-1,1,dz}];*)
(*][[1]]];*)
(*Export[FileNameJoin[{NotebookDirectory[],"refractiveIndex3d_tiltAngle_"<>ToString[DecimalForm[N@tiltAngle,{3,2}]]<>".dat"}],refractiveIndexData3d];*)


(* ::Input:: *)
(*xIndex=50;*)
(*plot=ListDensityPlot[Flatten[refractiveIndexData3d[[xIndex,All,All,{2,3,4}]],1],InterpolationOrder->0,FrameLabel->{"x","z"},PlotRange->All,PlotLegends->BarLegend[Automatic,LegendLabel->Placed["n",Right]](*,ColorFunctionScaling->False,ColorFunction->(ColorData["M10DefaultDensityGradient"][Rescale[#,{nOil,1.7}]]&)*),PlotRangePadding->None]*)
(*Histogram[Flatten[refractiveIndexData3d[[All,All,All,4]]],{nOil,1.02nLCe,0.01},Frame->True,FrameLabel->{"refractive index n","occurrence"},Prolog->{Gray,InfiniteLine[{{#,0},{#,1}}]&/@{nLCo,nLCe}}]*)


(* ::Input:: *)
(*dzz=10^6 rDropletDim dz; (** integral element **)*)
(*pathLengths3d=Table[Total[refractiveIndexData3d[[x,y,All,4]] dzz],{x,101},{y,101}];*)
(*Export[FileNameJoin[{NotebookDirectory[],"optical_pathLengths.dat"}],pathLengths3d];*)


(* ::Input:: *)
(*ListDensityPlot[pathLengths3d,InterpolationOrder->0,PlotRangePadding->None,PlotLegends->BarLegend[Automatic,LegendLabel->"optical pathlength (\[Mu]m)"],FrameLabel->{"x","y"}]*)
(**)
(*pathLengthCut=pathLengths3d[[xIndex,All]];*)
(*plot=ListLinePlot[{Rescale[Range[Length[pathLengthCut]],{1,Length[pathLengthCut]},{-1.0,1}],pathLengthCut}\[Transpose],Frame->True,FrameLabel->{"y (\!\(\*SubscriptBox[\(r\), \(droplet\)]\))","optical pathlength (\[Mu]m)"},FrameStyle->Directive[Black,14,AbsoluteThickness[2]],Axes->False]*)
(**)
(*pathLengthCut=pathLengths3d[[All,xIndex]];*)
(*plot=ListLinePlot[{Rescale[Range[Length[pathLengthCut]],{1,Length[pathLengthCut]},{-1.0,1}],pathLengthCut}\[Transpose],Frame->True,FrameLabel->{"x (\!\(\*SubscriptBox[\(r\), \(droplet\)]\))","optical pathlength (\[Mu]m)"},FrameStyle->Directive[Black,14,AbsoluteThickness[2]],Axes->False]*)


(* ::Subsection::Closed:: *)
(*save refractive index*)


(* ::Input:: *)
(*refIndexWater=1.333;*)
(*refIndexOil=1.433;*)
(*director[x_,y_,z_]=Normalize[{u[x,y,z],v[x,y,z],w[x,y,z]}/.solution];*)
(*refIndexLC[x_,y_,z_]=Norm[surfaceProjection[-{0,0,1},director[x,y,z]]];*)
(*lookupTablePositionRefractiveIndex=Flatten[Table[Flatten@{{x,y,z},If[Norm[{x,y,z}]<=1,If[z>0,refIndexLC[x,y,z],refIndexOil],refIndexWater]},{x,0.001,1,0.2},{y,-1,1,0.2},{z,-1,1,0.2}],2];*)
(*Export[FileNameJoin[{NotebookDirectory[],"positionRefractiveIndexLookupTable.dat"}],lookupTablePositionRefractiveIndex];*)
(*Export[FileNameJoin[{NotebookDirectory[],"positionRefractiveIndexLookupTable.mat"}],lookupTablePositionRefractiveIndex];*)


(* ::Input:: *)
(*dxData=Flatten[Table[Flatten@{{x,y,z},If[Norm[{x,y,z}]<=1,If[z>0,director[x,y,z][[1]],refIndexOil],refIndexWater]},{x,0.001,1,0.2},{y,-1,1,0.2},{z,-1,1,0.2}],2];*)


(* ::Input:: *)
(*ListSliceDensityPlot3D[dxData,"BackPlanes"(*"CenterPlanes"*),BoxRatios->{0.5,1,1}]*)


(* ::Input:: *)
(*ListSliceDensityPlot3D[lookupTablePositionRefractiveIndex,"CenterPlanes"]*)


(* ::Subsection::Closed:: *)
(*derivatives of director*)


(* ::Input:: *)
(*dxdirector[x_,y_,z_]=D[Evaluate[director[x,y,z]],x];*)


(* ::Input:: *)
(*dxdirector[xx_,yy_,zz_]=D[director[x,y,z],x]/.{x:>xx};*)
(*lookupTablePositionDirectorDx=Flatten[Table[{{x,y,z},If[Norm[{x,y,z}]<=1,dxdirector[x,y,z],{0,0,0}]},{x,0,1,0.2},{y,-1,1,0.2},{z,0.01,1,0.2}],2];*)
(*Export[FileNameJoin[{NotebookDirectory[],"dxdirectorLookupTable.dat"}],lookupTablePositionDirectorDx];*)
(**)
(*dydirector[x_,y_,z_]:=Derivative[0,1,0][Normalize[{u[x,y,z],v[x,y,z],w[x,y,z]}/.solution]];*)
(*lookupTablePositionDirectorDy=Flatten[Table[{{x,y,z},If[Norm[{x,y,z}]<=1,dydirector[x,y,z],{0,0,0}]},{x,0,1,0.2},{y,-1,1,0.2},{z,0.01,1,0.2}],2];*)
(*Export[FileNameJoin[{NotebookDirectory[],"dydirectorLookupTable.dat"}],lookupTablePositionDirectorDy];*)
(**)
(*dzdirector[x_,y_,z_]:=Derivative[0,0,1][Normalize[{u[x,y,z],v[x,y,z],w[x,y,z]}/.solution]];*)
(*lookupTablePositionDirectorDz=Flatten[Table[{{x,y,z},If[Norm[{x,y,z}]<=1,dzdirector[x,y,z],{0,0,0}]},{x,0,1,0.2},{y,-1,1,0.2},{z,0.01,1,0.2}],2];*)
(*Export[FileNameJoin[{NotebookDirectory[],"dzdirectorLookupTable.dat"}],lookupTablePositionDirectorDz];*)


(* ::Input:: *)
(*u[x,y,z]/.solution/.{x->0,y->0,z->0}*)
(*D[u[x,y,z]/.solution,x]/.{x->0,y->0,z->0}*)


(* ::Input:: *)
(*lookupTablePositionDirector=Flatten[Table[{{x,y,z},If[Norm[{x,y,z}]<=1,Normalize[{u[x,y,z],v[x,y,z],w[x,y,z]}/.solution],{0,0,0}]},{x,-1,1,0.2},{y,-1,1,0.2},{z,0.01,1,0.2}],2];*)
(*Export[FileNameJoin[{NotebookDirectory[],"directorLookupTable.mat"}],lookupTablePositionDirector];*)


(* ::Input:: *)
(*SliceDensityPlot3D[Evaluate[u[x,y,z]/.solution],{x,y,z}\[Element]halfBall,BoxRatios->{1,1,0.5},PlotLegends->Automatic]*)
(*SliceDensityPlot3D[Evaluate[w[x,y,z]/.solution],{x,y,z}\[Element]halfBall,BoxRatios->{1,1,0.5},PlotLegends->Automatic]*)


(* ::Subsection::Closed:: *)
(*Fresnel propagator*)


(* ::Subsubsection::Closed:: *)
(*Fresnel propagator test with square: 2d fourier transform*)


(* ::Input:: *)
(*eyTestSquare=Table[*)
(*If[Abs[i-50]<25\[And]Abs[j-50]<25,1,0]*)
(*,{i,nx},{j,ny}];*)
(*dx=0.1;*)
(*eyTestSin=Table[*)
(*Sin[i dx 2]Sin[j dx 2]*)
(*,{i,nx},{j,ny}];*)
(*ey=eyTestSin;*)


(* ::Input:: *)
(*ListDensityPlot[ey,InterpolationOrder->0,PlotRangePadding->None,FrameStyle->Directive[Black]]*)


(* ::Input:: *)
(*psd=Abs[Fourier[ey]]^2;*)
(*ListDensityPlot[Reverse[psd],PlotRange->All(*,ScalingFunctions->"Log"*),PlotRangePadding->None]*)


(* ::Input:: *)
(*Clear[\[Lambda]]*)


(* ::Input:: *)
(*Exp[I 2\[Pi] z/\[Lambda]]/(I \[Lambda] z)FourierTransform[Exp[I \[Pi] x^2/(\[Lambda] z)],x,kx]FourierTransform[Exp[I \[Pi] y^2/(\[Lambda] z)],y,ky]//Simplify*)


(* ::Input:: *)
(*E^(I k  z-(I (kx^2+ky^2) z \[Lambda])/(4 \[Pi]))/(2 \[Pi])*)


(* ::Input:: *)
(*d=d*(-1)^Table[i+j,{i,nRow},{j,nCol}];  (*center,i.e.fftshift like*)*)


(* ::Input:: *)
(*Exp[I 2\[Pi] z/\[Lambda]]/(I \[Lambda] z)FourierTransform[Exp[I \[Pi] (x^2+y^2)/(\[Lambda] z)],{x,y},{kx,ky}]*)


(* ::Input:: *)
(*fresnelPropagatorKernelFourier[z_,dkx_,dky_,\[Lambda]_]:=Table[E^((2 I \[Pi] z)/\[Lambda]-(I (kx^2+ky^2) z \[Lambda])/(4 \[Pi]))/(2 \[Pi]),{i,nx},{j,ny}];*)


(* ::Input:: *)
(*{dx,dy}=10^9 2rDropletDim/nx{1,1};  (** nm **)*)
(*p0=ListDensityPlot[Abs[ey]^2,InterpolationOrder->0,PerformanceGoal->"Quality",ColorFunction->GrayLevel,PlotRangePadding->None];*)
(*Manipulate[*)
(*zz=a (*1.2*)(*1000*) 10^9 \[Lambda]; (** nm **)*)
(*kernel=fresnelPropagatorKernel[zz,dx,dy,10^9 \[Lambda]];*)
(*eySmooth=GaussianFilter[ey,{5,0}];*)
(*propagatedEy=ListConvolve[kernel,eySmooth,{1,1}];  (** periodic **)*)
(*(*propagatedEy=ListConvolve[kernel,eySmooth,{1,-1},{0,0}];*) (** zero padding **)*)
(**)
(*p1=ListDensityPlot[Abs[eySmooth]^2,InterpolationOrder->0,PerformanceGoal->"Quality",ColorFunction->GrayLevel,PlotRangePadding->None];*)
(*p2=ListDensityPlot[Re[kernel],InterpolationOrder->0,PerformanceGoal->"Quality",ColorFunction->GrayLevel,PlotRange->All,PlotRangePadding->None];*)
(*p3=ListDensityPlot[Abs[propagatedEy(*\[LeftDoubleBracket]51;;150,51;;150\[RightDoubleBracket]*)]^2,InterpolationOrder->0,(*PlotLegends->Automatic,*)PlotRange->All,PlotRangePadding->None,FrameStyle->Directive[Black],(*FrameLabel->{"x","y"},*)PerformanceGoal->"Quality",ColorFunction->GrayLevel];*)
(*Grid[{{(*p0,p1,*)p2,p3}}]*)
(*,{{a,1,"multiplication factor"},1,100,0.1,Appearance->"Open"},TrackedSymbols:>{a}]*)


(* ::Subsubsection::Closed:: *)
(*Fresnel propagator test with square: 2d real space convolution*)


(* ::Input:: *)
(*eyTestSquare=Table[*)
(*If[Abs[i-50]<25\[And]Abs[j-50]<25,1,0]*)
(*,{i,nx},{j,ny}];*)


(* ::Input:: *)
(*ListDensityPlot[eyTestSquare,InterpolationOrder->0]*)


(* ::Input:: *)
(*data=N@Table[Exp[-((i -50)/10)^2],{i,100}];*)
(*ListPlot[data,PlotRange->All]*)
(*convData=ListConvolve[{1}(*ConstantArray[1,4(*Length[data]*)]*),data,1];*)
(*ListPlot[convData,PlotRange->All]*)


(* ::Input:: *)
(*convData2d=ListConvolve[{{1,0},{0,0}},eySmooth,{1,1}];*)


(* ::Input:: *)
(*ListDensityPlot[convData2d]*)


(* ::Input:: *)
(*ey=eyTestSquare;*)


(* ::Input:: *)
(*{dx,dy}=10^9 2rDropletDim/nx{1,1};  (** nm **)*)
(*p0=ListDensityPlot[Abs[ey]^2,InterpolationOrder->0,PerformanceGoal->"Quality",ColorFunction->GrayLevel,PlotRangePadding->None];*)
(*Manipulate[*)
(*zz=a (*1.2*)(*1000*) 10^9 \[Lambda]; (** nm **)*)
(*kernel=fresnelPropagatorKernel[zz,dx,dy,10^9 \[Lambda]];*)
(*eySmooth=GaussianFilter[ey,{5,0}];*)
(*propagatedEy=ListConvolve[kernel,eySmooth,{1,1}];  (** periodic **)*)
(*(*propagatedEy=ListConvolve[kernel,eySmooth,{1,-1},{0,0}];*) (** zero padding **)*)
(**)
(*p1=ListDensityPlot[Abs[eySmooth]^2,InterpolationOrder->0,PerformanceGoal->"Quality",ColorFunction->GrayLevel,PlotRangePadding->None];*)
(*p2=ListDensityPlot[Re[kernel],InterpolationOrder->0,PerformanceGoal->"Quality",ColorFunction->GrayLevel,PlotRange->All,PlotRangePadding->None];*)
(*p3=ListDensityPlot[Abs[propagatedEy(*\[LeftDoubleBracket]51;;150,51;;150\[RightDoubleBracket]*)]^2,InterpolationOrder->0,(*PlotLegends->Automatic,*)PlotRange->All,PlotRangePadding->None,FrameStyle->Directive[Black],(*FrameLabel->{"x","y"},*)PerformanceGoal->"Quality",ColorFunction->GrayLevel];*)
(*Grid[{{(*p0,p1,*)p2,p3}}]*)
(*,{{a,1,"multiplication factor"},1,100,0.1,Appearance->"Open"},TrackedSymbols:>{a}]*)


(* ::Subsubsection::Closed:: *)
(*Fresnel propagator: 2d real space convolution*)


(* ::Input:: *)
(*Import["eyRGB.mat",ey]*)


(* ::Input:: *)
(*ex=afterDropletEField[[All,All,1]];*)
(*ey=afterDropletEField[[All,All,2]];*)


(* ::Input:: *)
(*ListDensityPlot[Abs@ey]*)


(* ::Input:: *)
(*Fourier*)


fresnelPropagatorKernel[z_,dx_,dy_,\[Lambda]_]:=Table[
Exp[I 2\[Pi] z/\[Lambda]]/(I \[Lambda] z) Exp[I \[Pi] (((i-0.5nx)dx)^2+(((j -0.5ny)dy))^2)/(\[Lambda] z)]
,{i,nx},{j,ny}];


(* ::Input:: *)
(*{dx,dy}=10^9 2rDropletDim/nx{1,1} (** nm **)*)
(*Manipulate[*)
(*(*zz=a (*1.2*)(*1000*) 10^9 \[Lambda]; (** nm **)*)*)
(**)
(*kernel=fresnelPropagatorKernel[zz,dx,dy,10^9 \[Lambda]];*)
(*(*eySmooth=GaussianFilter[ey,{5,0}];*)*)
(*eySmooth=ey;*)
(*propagatedEy=ListConvolve[kernel,eySmooth,{1,-1}];  (** periodic **)*)
(*(*propagatedEy=ListConvolve[kernel,eySmooth,{1,-1},{0,0}];*) (** zero padding **)*)
(**)
(*(*p0=ListDensityPlot[Abs[ey]^2,InterpolationOrder->0,PerformanceGoal->"Quality",ColorFunction->GrayLevel,PlotRangePadding->None];*)
(*p1=ListDensityPlot[Abs[eySmooth]^2,InterpolationOrder->0,PerformanceGoal->"Quality",ColorFunction->GrayLevel,PlotRangePadding->None];*)*)
(*p2=ListDensityPlot[Re[kernel],InterpolationOrder->0,PerformanceGoal->"Quality",ColorFunction->GrayLevel,PlotRange->All,PlotRangePadding->None];*)
(*p3=ListDensityPlot[Abs[propagatedEy[[51;;150,51;;150]]]^2,InterpolationOrder->0,(*PlotLegends->Automatic,*)PlotRange->All,PlotRangePadding->None,FrameStyle->Directive[Black],(*FrameLabel->{"x","y"},*)PerformanceGoal->"Quality",ColorFunction->GrayLevel];*)
(*Grid[{{(*p0,*)(*p1,*)p2,p3}}]*)
(*(*,{{a,1,"multiplication factor"},1,100,0.1,Appearance->"Open"}*)*)
(*,{{zz,1,"z distance in nm"},1,10000,0.1,Appearance->"Open"}*)
(*,TrackedSymbols:>{zz}]*)


(* ::Chapter:: *)
(*join RGB channels together*)


(* ::Input:: *)
(*images=Import[FileNameJoin[{ParentDirectory[NotebookDirectory[]],"Droplet optics",#,"zDistance_0.0","all_intensity_images_overview.png"}]]&/@names;*)


(* ::Input:: *)
(*grayImages[[1]]==grayImages[[2]]*)


(* ::Input:: *)
(*intensityRGBdata//Dimensions*)


(* ::Input:: *)
(*intensityRGBdata=ImageData[ColorConvert[#,"Grayscale"]]&/@images;*)


(* ::Input:: *)
(*Image[intensityRGBdata[[All,1;;1000,1;;1000]],"Byte",Interleaving->False]*)


(* ::Input:: *)
(*names={*)
(*"prop_test_dropletTilt_thetaPhi_analyticDirector_wvl_630_DnLC_0.17",*)
(*"prop_test_dropletTilt_thetaPhi_analyticDirector_wvl_550_DnLC_0.18",*)
(*"prop_test_dropletTilt_thetaPhi_analyticDirector_wvl_470_DnLC_0.20"*)
(*};*)
(*thetaPhiIntensitiesRGB=Import[FileNameJoin[{ParentDirectory[NotebookDirectory[]],"Droplet optics",#,"zDistance_0.0","integratedIntensity.mat"}]][[1]]&/@names;*)


(* ::Input:: *)
(*thetaPhiIntensitiesRGB//Dimensions*)
(*meanIntensities=Mean[thetaPhiIntensitiesRGB];*)
(*meanIntensities=Rescale[meanIntensities];*)
(*Export[FileNameJoin[{ParentDirectory[NotebookDirectory[]],"Droplet optics","RGB","normIntegratedIntensity.mat"}],meanIntensities]*)


(* ::Input:: *)
(**)


(* ::Input:: *)
(*Plot3D[2000+1500Sin[2\[Theta]](1+Sin[2\[Phi]]),{\[Theta],0,\[Pi]/2},{\[Phi],0,\[Pi]/4},ColorFunction->"BlueGreenYellow"]*)
(*ListPlot3D[meanIntensities\[Transpose],DataRange->{{0,90},{0,45}},ColorFunction->"BlueGreenYellow",PlotRangePadding->None]*)


(* ::Input:: *)
(*ListContourPlot[meanIntensities\[Transpose],DataRange->{{0,90},{0,45}},ColorFunction->"BlueGreenYellow",PlotRangePadding->None]*)
(*ListContourPlot[thetaPhiIntensitiesRGB[[1]]\[Transpose],DataRange->{{0,90},{0,45}},ColorFunction->"BlueGreenYellow",PlotRangePadding->None]*)
(*ListContourPlot[thetaPhiIntensitiesRGB[[2]]\[Transpose],DataRange->{{0,90},{0,45}},ColorFunction->"BlueGreenYellow",PlotRangePadding->None]*)
(*ListContourPlot[thetaPhiIntensitiesRGB[[3]]\[Transpose],DataRange->{{0,90},{0,45}},ColorFunction->"BlueGreenYellow",PlotRangePadding->None]*)


(* ::Chapter:: *)
(*Figure S1c*)


(* ::Input:: *)
(*name="RGB";*)
(*updateIntensityLookupMap[name];*)


(* ::Input:: *)
(*(** \[Theta]=radial coordinate, \[Phi]=angle coordinate **)*)
(*data=Flatten[Table[{\[Theta] Cos[\[Phi]],\[Theta] Sin[\[Phi]],thetaPhiIntensityInterpolation[-Abs[\[Theta]-\[Pi]/2.0]+\[Pi]/2.0,-Abs[Mod[\[Phi],\[Pi]/2.0]-\[Pi]/4.0]+\[Pi]/4.0]},{\[Theta],Subdivide[0,\[Pi]/2,100-1]},{\[Phi],Most[Subdivide[0,2\[Pi],100-1]]}],1];*)
(*plot1c=ListDensityPlot[data,PlotRange->All,ColorFunction->GrayLevel,Epilog->{Darker@Red,*)
(*Arrow[{{0,0},{\[Pi]/2,0}}],Text[Style["\[Theta]",18],Scaled[{0.92,0.5}]],Arrow@ResourceFunction["SplineCircle"][{0,0},1.1\[Pi]/2,{1,0},{5Degree,90 Degree}],Text[Style["\[Phi]",18],Scaled[{0.83,0.83}]]},Frame->False,PlotRangePadding->Scaled[0.1]]*)
(**)
(*Export[FileNameJoin[{"C:\\Users\\Jan\\Dropbox (MIT)\\Bacteria Droplet Manuscript\\Figures\\Figures current\\SI figures\\intensity","SI_figure_1c.png"}],plot1c]*)


(* ::Chapter::Closed:: *)
(*old*)


(* ::Subsubsection::Closed:: *)
(*compute LC alignment in upper half of droplet for raytracer simulation*)


(* ::Text:: *)
(*W. Wang, L. Zhang, P. Zhang, Modelling and computation of liquid crystals. Acta Numerica. 30, 765\[Dash]851 (2021)*)


(* ::Input:: *)
(*Tr[Grad[{u[x,y,z],v[x,y,z],w[x,y,z]},{x,y,z}]\[Transpose]Grad[{u[x,y,z],v[x,y,z],w[x,y,z]},{x,y,z}]]*)


(* ::Input:: *)
(*Clear[u,x,y,z];*)
(*halfBall=ImplicitRegion[x^2+y^2+z^2<=1\[And]z>=0,{x,y,z}];*)
(*bcSphere=(IdentityMatrix[3]-TensorProduct[{x,y,z},{x,y,z}]) . {0,0,1};*)
(*solution=NDSolve[{*)
(*Laplacian[u[x,y,z],{x,y,z}]==-(\!\(\*SuperscriptBox[\(w\), *)
(*TagBox[*)
(*RowBox[{"(", *)
(*RowBox[{"0", ",", "0", ",", "1"}], ")"}],*)
(*Derivative],*)
(*MultilineFunction->None]\)[x,y,z]^2+\!\(\*SuperscriptBox[\(v\), *)
(*TagBox[*)
(*RowBox[{"(", *)
(*RowBox[{"0", ",", "1", ",", "0"}], ")"}],*)
(*Derivative],*)
(*MultilineFunction->None]\)[x,y,z]^2+\!\(\*SuperscriptBox[\(u\), *)
(*TagBox[*)
(*RowBox[{"(", *)
(*RowBox[{"1", ",", "0", ",", "0"}], ")"}],*)
(*Derivative],*)
(*MultilineFunction->None]\)[x,y,z]^2)u[x,y,z],*)
(*Laplacian[v[x,y,z],{x,y,z}]==-(\!\(\*SuperscriptBox[\(w\), *)
(*TagBox[*)
(*RowBox[{"(", *)
(*RowBox[{"0", ",", "0", ",", "1"}], ")"}],*)
(*Derivative],*)
(*MultilineFunction->None]\)[x,y,z]^2+\!\(\*SuperscriptBox[\(v\), *)
(*TagBox[*)
(*RowBox[{"(", *)
(*RowBox[{"0", ",", "1", ",", "0"}], ")"}],*)
(*Derivative],*)
(*MultilineFunction->None]\)[x,y,z]^2+\!\(\*SuperscriptBox[\(u\), *)
(*TagBox[*)
(*RowBox[{"(", *)
(*RowBox[{"1", ",", "0", ",", "0"}], ")"}],*)
(*Derivative],*)
(*MultilineFunction->None]\)[x,y,z]^2)v[x,y,z],*)
(*Laplacian[w[x,y,z],{x,y,z}]==-(\!\(\*SuperscriptBox[\(w\), *)
(*TagBox[*)
(*RowBox[{"(", *)
(*RowBox[{"0", ",", "0", ",", "1"}], ")"}],*)
(*Derivative],*)
(*MultilineFunction->None]\)[x,y,z]^2+\!\(\*SuperscriptBox[\(v\), *)
(*TagBox[*)
(*RowBox[{"(", *)
(*RowBox[{"0", ",", "1", ",", "0"}], ")"}],*)
(*Derivative],*)
(*MultilineFunction->None]\)[x,y,z]^2+\!\(\*SuperscriptBox[\(u\), *)
(*TagBox[*)
(*RowBox[{"(", *)
(*RowBox[{"1", ",", "0", ",", "0"}], ")"}],*)
(*Derivative],*)
(*MultilineFunction->None]\)[x,y,z]^2)w[x,y,z],*)
(*DirichletCondition[u[x,y,z]==bcSphere[[1]],x^2+y^2+z^2==1],*)
(*DirichletCondition[v[x,y,z]==bcSphere[[2]],x^2+y^2+z^2==1],*)
(*DirichletCondition[w[x,y,z]==bcSphere[[3]],x^2+y^2+z^2==1],*)
(*DirichletCondition[u[x,y,z]==0,z==0],*)
(*DirichletCondition[v[x,y,z]==0,z==0],*)
(*DirichletCondition[w[x,y,z]==1,z==0]*)
(*},{u[x,y,z],v[x,y,z],w[x,y,z]},{x,y,z}\[Element]halfBall][[1]];*)


(* ::Input:: *)
(*director[x_,y_,z_]=Normalize[{u[x,y,z],v[x,y,z],w[x,y,z]}/.solution];*)
(*lookupTablePositionDirector=Flatten[Table[{{x,y,z},If[Norm[{x,y,z}]<=1,director[x,y,z],{0,0,0}]},{x,0,1,0.2},{y,-1,1,0.2},{z,0.01,1,0.2}],2];*)
(*Export[FileNameJoin[{NotebookDirectory[],"directorLookupTable.dat"}],lookupTablePositionDirector];*)


(* ::Input:: *)
(*refIndexWater=1.333;*)
(*director[x_,y_,z_]=Normalize[{u[x,y,z],v[x,y,z],w[x,y,z]}/.solution];*)
(*lookupTablePositionRefractiveIndex=Flatten[Table[Flatten@{{x,y,z},If[Norm[{x,y,z}]<=1,Norm[surfaceProjection[-{0,0,1},director[x,y,z]]],refIndexWater]},{x,0,1,0.2},{y,-1,1,0.2},{z,0.01,1,0.2}],2];*)
(*Export[FileNameJoin[{NotebookDirectory[],"positionRefractiveIndexLookupTable.dat"}],lookupTablePositionRefractiveIndex];*)
(*Export[FileNameJoin[{NotebookDirectory[],"positionRefractiveIndexLookupTable.mat"}],lookupTablePositionRefractiveIndex];*)


(* ::Input:: *)
(*director[x,y,z]*)


(* ::Input:: *)
(*dxdirector[x_,y_,z_]=D[Evaluate[director[x,y,z]],x];*)


(* ::Input:: *)
(*director[0,0,0]*)


(* ::Input:: *)
(*dxdirector[0.1,0.1,0.1]*)


(* ::Input:: *)
(*dxdirector[xx_,yy_,zz_]=D[director[x,y,z],x]/.{x:>xx};*)
(*lookupTablePositionDirectorDx=Flatten[Table[{{x,y,z},If[Norm[{x,y,z}]<=1,dxdirector[x,y,z],{0,0,0}]},{x,0,1,0.2},{y,-1,1,0.2},{z,0.01,1,0.2}],2];*)
(*Export[FileNameJoin[{NotebookDirectory[],"dxdirectorLookupTable.dat"}],lookupTablePositionDirectorDx];*)
(**)
(*dydirector[x_,y_,z_]:=Derivative[0,1,0][Normalize[{u[x,y,z],v[x,y,z],w[x,y,z]}/.solution]];*)
(*lookupTablePositionDirectorDy=Flatten[Table[{{x,y,z},If[Norm[{x,y,z}]<=1,dydirector[x,y,z],{0,0,0}]},{x,0,1,0.2},{y,-1,1,0.2},{z,0.01,1,0.2}],2];*)
(*Export[FileNameJoin[{NotebookDirectory[],"dydirectorLookupTable.dat"}],lookupTablePositionDirectorDy];*)
(**)
(*dzdirector[x_,y_,z_]:=Derivative[0,0,1][Normalize[{u[x,y,z],v[x,y,z],w[x,y,z]}/.solution]];*)
(*lookupTablePositionDirectorDz=Flatten[Table[{{x,y,z},If[Norm[{x,y,z}]<=1,dzdirector[x,y,z],{0,0,0}]},{x,0,1,0.2},{y,-1,1,0.2},{z,0.01,1,0.2}],2];*)
(*Export[FileNameJoin[{NotebookDirectory[],"dzdirectorLookupTable.dat"}],lookupTablePositionDirectorDz];*)


(* ::Input:: *)
(*u[x,y,z]/.solution/.{x->0,y->0,z->0}*)
(*D[u[x,y,z]/.solution,x]/.{x->0,y->0,z->0}*)


(* ::Input:: *)
(*lookupTablePositionDirector=Flatten[Table[{{x,y,z},If[Norm[{x,y,z}]<=1,Normalize[{u[x,y,z],v[x,y,z],w[x,y,z]}/.solution],{0,0,0}]},{x,-1,1,0.2},{y,-1,1,0.2},{z,0.01,1,0.2}],2];*)
(*Export[FileNameJoin[{NotebookDirectory[],"directorLookupTable.mat"}],lookupTablePositionDirector];*)


(* ::Input:: *)
(*points=Flatten[*)
(*Table[Table[*)
(*angle=ArcSin[z];*)
(*N[{#[[1]],#[[2]],Sin[angle]}&/@(0.95r Cos[angle]CirclePoints[Round[15r Cos[angle]+5]])],{r,Subdivide[0.0,0.99,Round[10(1-z)]]}],{z,Subdivide[0.1,0.95,3]}],2];*)
(*vectors=Table[*)
(*{x,y,z}=p;*)
(*orientation=Normalize[{u[x,y,z],v[x,y,z],w[x,y,z]}/.solution];*)
(*({x,y,z}+# 0.05orientation)&/@{-1,1}*)
(*,{p,points}];*)


(* ::Input:: *)
(*plot=Graphics3D[{Evaluate[CapsuleShape[#,0.03]&/@vectors],Opacity[0.2],Sphere[]},PlotRange->{All,All,{0,All}},Boxed->False,ViewVertical->{-0.11293301187983383`,-0.07449068114514881`,0.9908063752571847`},ViewPoint->{-2.4149686489954596`,-1.5612429571636024`,1.7833807369925514`}]*)
(*Export[FileNameJoin[{NotebookDirectory[],"directorDistribution.png"}],plot]*)


(* ::Input:: *)
(*SliceDensityPlot3D[Evaluate[u[x,y,z]/.solution],{x,y,z}\[Element]halfBall,BoxRatios->{1,1,0.5},PlotLegends->Automatic]*)
(*SliceDensityPlot3D[Evaluate[w[x,y,z]/.solution],{x,y,z}\[Element]halfBall,BoxRatios->{1,1,0.5},PlotLegends->Automatic]*)


(* ::Subsection::Closed:: *)
(*malteser cross as a function of tilt angles \[Theta], \[Phi] with variable polarization*)


(* ::Subsubsection::Closed:: *)
(*main*)


(* ::Input:: *)
(*{nx,ny,nz}={100,100,100};*)
(*xVals=Subdivide[-1.0,1.0,nx-1];*)
(*yVals=Subdivide[-1.0,1.0,ny-1];*)
(*zVals=Subdivide[-1.0,1.0,nz-1];*)
(*dz=(zVals[[2]]-zVals[[1]]);*)
(*{fcAngle,waterAngle}={-20(*10*),-10};*)
(**)
(*d=rDropletDim dz;*)
(*\[Delta]LC[nLCe_,nLCo_]:=2\[Pi]/\[Lambda] d (nLCe-nLCo);*)
(*\[Delta]o[nLCo_]:=2\[Pi]/\[Lambda] d nLCo;*)
(*(*Inverse[RotationMatrix[\[Theta]]].{{1,0},{0,Exp[\[ImaginaryI] \[Delta]]}}.RotationMatrix[\[Theta]]//Simplify*)*)
(*jonesMatrixLC[\[Theta]_,nLCe_]:=Exp[I \[Delta]o[nLCo]]{{Cos[\[Theta]]^2+Sin[\[Theta]]^2 Exp[I \[Delta]LC[nLCe,nLCo]],Sin[\[Theta]]Cos[\[Theta]](Exp[I \[Delta]LC[nLCe,nLCo]]-1)},{Sin[\[Theta]]Cos[\[Theta]](Exp[I \[Delta]LC[nLCe,nLCo]]-1),Sin[\[Theta]]^2+Cos[\[Theta]]^2 Exp[I \[Delta]LC[nLCe,nLCo]]}};*)
(*jonesMatrixFC=Exp[I \[Delta]o[nOil]]IdentityMatrix[2];*)
(*jonesMatrixFCdroplet[d_]:=Exp[I 2\[Pi]/\[Lambda] d nOil]IdentityMatrix[2];*)
(*jonesMatrixVerticalPolarizer={{0,0},{0,1}};*)
(**)
(*tiltAnglesTheta=N[Subdivide[0,\[Pi]/2,20]];*)
(*tiltAnglesPhi=N[Subdivide[0,\[Pi]/4,10]];*)
(*tiltAnglesTheta={0.0};*)
(*tiltAnglesPhi={0.0};*)
(*tiltDirector[\[Theta]_,\[Phi]_]:={Cos[\[Phi]]Sin[\[Theta]],Sin[\[Phi]]Sin[\[Theta]],Cos[\[Theta]]};*)
(*saveStateQ=False;*)
(*Print[AbsoluteTiming[*)
(*intensityTiltData=Table[*)
(*rotaMatrix=RotationMatrix[{{0,0,1},tiltDirector[tiltAngleTheta,tiltAnglePhi]}];*)
(*rotaMatrixInv=RotationMatrix[{tiltDirector[tiltAngleTheta,tiltAnglePhi],{0,0,1}}];*)
(**)
(*(** LC director projected on xy plane angle = marker for liquid type (LC,FC,water) **)*)
(*(** calculate extraordinary refractive index Subscript[n, e]'' for LC director != wave director **)*)
(*zdirectorAngles3d=Table[*)
(*position={x,y,z};*)
(*directorAngle=If[Norm[position]<=0.99,*)
(*If[position . tiltDirector[tiltAngleTheta,tiltAnglePhi]>0.0,*)
(*rotatedPosition=rotaMatrixInv . position;*)
(*rotatedDirector=rotaMatrix . director[rotatedPosition];*)
(*lcAngle=ArcTan[rotatedDirector[[1]],rotatedDirector[[2]]];*)
(*lcAngle,fcAngle],waterAngle];*)
(*{z,directorAngle}*)
(*,{x,xVals},{y,yVals},{z,zVals}];*)
(**)
(*(*nLC=calcLCrefIndex[rotatedDirector,polarizationVector];*)*)
(*(** calculate polarization state of light in LC part + FC part **)*)
(*dropletPolarizationState=Table[*)
(*zdirectorAngles=Reverse@zdirectorAngles3d[[xIndex,yIndex]];*)
(*filteredzAngles=If[#[[2]]!=waterAngle,#,Nothing]&/@zdirectorAngles;*)
(*afterVerticalPolarizerState=If[filteredzAngles!={},*)
(*(** ray intersection with the droplet **)*)
(*polarizationState=polarizationStateHorizontal;*)
(*lczAngles=If[#[[2]]!=fcAngle,#,Nothing]&/@filteredzAngles;*)
(*If[lczAngles!={},*)
(*Do[*)
(*{z,lcAngle}=lczAngle;*)
(*position={xVals[[xIndex]],yVals[[yIndex]],z};*)
(*rotatedPosition=rotaMatrixInv . position;*)
(*rotatedDirector=rotaMatrix . director[rotatedPosition];*)
(*nLCeMod=calcLCrefIndex[rotatedDirector,{#[[1]],#[[2]],0.0}&@polarizationState];*)
(*polarizationState=jonesMatrixLC[lcAngle,nLCeMod] . polarizationState;*)
(*,{lczAngle,lczAngles}];*)
(*];*)
(*afterLCState=polarizationState;*)
(*fcBoxCount=Count[filteredzAngles[[All,2]],fcAngle];*)
(*afterDropletState=jonesMatrixFCdroplet[dz fcBoxCount] . afterLCState;*)
(*jonesMatrixVerticalPolarizer . afterDropletState,*)
(*(** no ray intersection with the droplet **)*)
(*jonesMatrixVerticalPolarizer . polarizationStateHorizontal*)
(*];*)
(*afterVerticalPolarizerState*)
(*,{xIndex,Dimensions[zdirectorAngles3d][[1]]},{yIndex,Dimensions[zdirectorAngles3d][[2]]}];*)
(*intensity=Map[Abs[#[[1]]]^2+Abs[#[[2]]]^2&,dropletPolarizationState,{2}];*)
(**)
(*(** output **)*)
(*{tiltAngleTheta,tiltAnglePhi,Mean[Flatten[intensity]],dropletPolarizationState,intensity}*)
(*,{tiltAngleTheta,tiltAnglesTheta}*)
(*,{tiltAnglePhi,tiltAnglesPhi}];*)
(*][[1]]];*)
(**)
(*(** save polarization state **)*)
(*If[saveStateQ,*)
(*Do[*)
(*finalPolarizationState=Map[ReIm,intensityTiltData[[i,j,4]]];*)
(*Export[FileNameJoin[{$HomeDirectory,"Desktop","afterDropletPolarizationState_i_"<>ToString[i]<>"_j_"<>ToString[j]<>".mat"}],finalPolarizationState];*)
(*,{i,Length[tiltAnglesTheta]},{j,Length[tiltAnglesTheta]}]*)
(*];*)


(* ::Input:: *)
(*avgnLC=Table[*)
(*Mean[directorAnglesRefIndices3d[[xIndex,yIndex,All,2]]]*)
(*,{xIndex,Dimensions[directorAnglesRefIndices3d][[1]]}*)
(*,{yIndex,Dimensions[directorAnglesRefIndices3d][[2]]}];*)
(*ListDensityPlot[avgnLC,InterpolationOrder->0]*)


(* ::Input:: *)
(*ListContourPlot[avgnLC,Contours->{nz nWater,1400,1410(*nz(0.5 nLCo+0.5nLCe)*)}]*)


(* ::Input:: *)
(*plot=Show[*)
(*Plot[0.35Cos[\[Theta]]^2,{\[Theta],0,\[Pi]/2},PlotStyle->LightGray,PlotLegends->{"\!\(\*SuperscriptBox[\(cos\), \(2\)]\)\[Theta]"}],*)
(*ListLinePlot[intensityTiltData[[All,{1,2}]],Mesh->All,PlotLegends->{"Jones matrix"}]*)
(*,Frame->True,FrameLabel->{"tilt angle \[Theta]","intensity"},PlotLabel->Style["\[Phi] = "<>ToString[Round[tiltAnglePhi 180/\[Pi]]]<>"\[Degree]",Black,20],FrameStyle->Directive[Black,18,AbsoluteThickness[2]],AspectRatio->1*)
(*]*)
(*Export[FileNameJoin[{NotebookDirectory[],"intensity_over_tilt_phi_"<>ToString[IntegerString[tiltAnglePhi 180/\[Pi]]]<>".png"}],plot]*)


(* ::Input:: *)
(*i=-1;*)
(*ListDensityPlot[intensityTiltData[[1,1,5]],InterpolationOrder->0,PlotRangePadding->None,PlotRange->All,PlotLegends->BarLegend[Automatic,LegendLabel->"intensity"],FrameLabel->{"x","y"},PlotLabel->"\[Theta] = "<>ToString[DecimalForm[intensityTiltData[[i,1,1]]180/\[Pi],{4,1}]]<>" \[Degree]"*)
(*,ColorFunction->GrayLevel*)
(*(*,PerformanceGoal->"Speed"*)]*)


(* ::Input:: *)
(*counter=0;*)
(*Do[*)
(*label="\[Theta] = "<>ToString[DecimalForm[intensityTiltData[[i,1]]180/\[Pi],{4,1}]]<>"\[Degree], \[Phi] = "<>ToString[DecimalForm[tiltAnglePhi 180/\[Pi],{4,1}]]<>"\[Degree]";*)
(*plotHighRes=ListDensityPlot[intensityTiltData[[i,3]],InterpolationOrder->0,PlotRangePadding->None,PlotRange->All,PlotLegends->BarLegend[Automatic,LegendLabel->"intensity"],FrameLabel->{"x","y"},ColorFunction->GrayLevel,ImageSize->{Automatic,400}];*)
(*plotLowRes=ListDensityPlot[intensityTiltData[[i,3]],InterpolationOrder->0,PlotRangePadding->None,PlotRange->All,PlotLegends->BarLegend[Automatic,LegendLabel->"intensity"],FrameLabel->{"x","y"},ColorFunction->GrayLevel,ImageSize->{Automatic,400},PerformanceGoal->"Speed"];*)
(*plot=Grid[{{label,SpanFromLeft},{plotHighRes,plotLowRes}}];*)
(*Export[FileNameJoin[{NotebookDirectory[],"mma_polarizationMovie_"<>ToString[tiltAnglePhi 180/\[Pi]]<>"_deg","p_"<>IntegerString[++counter,10,IntegerLength[Length@intensityTiltData]+1]<>".png"}],plot];*)
(*,{i,Length[intensityTiltData]}];*)


(* ::Input:: *)
(*Manipulate[*)
(*ListDensityPlot[intensityTiltData[[i,3]],InterpolationOrder->0,PlotRangePadding->None,PlotRange->All,PlotLegends->BarLegend[Automatic,LegendLabel->"intensity"],FrameLabel->{"x","y"},PlotLabel->"\[Theta] = "<>ToString@intensityTiltData[[i,1]],ColorFunction->GrayLevel]*)
(*,{i,1,21,1}]*)
