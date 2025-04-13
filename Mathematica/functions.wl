(* ::Package:: *)

(* ::Title:: *)
(*Functions for Running and Tumbling on a Sphere*)


(* ::Chapter:: *)
(*init*)


(* ::Section::Closed:: *)
(*compiler settings*)


(* ::Text:: *)
(*Note: You need to have a C++ compiler installed on your system for this to work.*)
(*Examples: Visual Studio Community, mingw, ...*)


Needs["CCompilerDriver`GenericCCompiler`"];


If[($UserName=="Jan")\[And]($OperatingSystem=="Windows")\[And]($MachineName=="lepton"\[Or]$MachineName=="myon"),
	$CCompiler={
		"Compiler"->GenericCCompiler,
		"CompilerInstallation"->"C:\\msys64\\mingw64\\bin",
		"CompilerName"->"gcc.exe"
	};
];


(* ::Chapter:: *)
(*Utility functions*)


(* ::Section:: *)
(*general analysis*)


(* ::Subsection::Closed:: *)
(*angle2d*)


angle2d[v_,w_]:=ArcTan[v . w,v[[1]]w[[2]]-w[[1]]v[[2]]]


(* ::Subsection::Closed:: *)
(*calcMSD*)


(* ::Input:: *)
(*{msd,msdExponent}=calcMSD[trajectory,pthout];*)


(* ::Subsubsection:: *)
(*code*)


msdC=Compile[{{trajectory,_Real,2},{lag,_Integer,0}},
	Total[Mean/@Transpose[Differences[trajectory,1,lag]^2]]
,CompilationTarget->"C",RuntimeOptions->"Speed"];


Clear[calcMSD];
calcMSD[trajectory_]:=Module[{msd,msdData,y0,msdExponent(*,msdPlot,msdPlotLogLog*)},

(** (1) distance over varying time intervals to start point **)
msdData=Table[{lag,msdC[trajectory,lag]},{lag,1,Length[trajectory]-2}];

(** 2 plot things and get diffusion coefficient **)
{y0,msdExponent}=LinearModelFit[Log10[msdData],lag,lag]["ParameterTableEntries"][[All,1]];
msdPlot=ListPlot[msdData,FrameLabel->{"time interval \[CapitalDelta]t (s)","MSD(\[CapitalDelta]t)"},Frame->True,PlotRangePadding->None,FrameStyle->Directive[Black,20,FontFamily->"Helvetica",AbsoluteThickness[2]],ImageSize->500];
msdPlotLogLog=Show[
	ListLogLogPlot[msdData,FrameLabel->{"log(\[CapitalDelta]t)","log(MSD)"},Frame->True,PlotRangePadding->None,FrameStyle->Directive[Black,20,FontFamily->"Helvetica",AbsoluteThickness[2]],ImageSize->500,PlotLegends->"exponent: "<>ToString[msdExponent]],
	LogLogPlot[10^(msdExponent Log10[t]+y0),{t,0.00001,Length[trajectory]-2},PlotStyle->Red]
];
(*Export[FileNameJoin[{pthout,"msdplot.png"}],msdPlot];
Export[FileNameJoin[{pthout,"msdplotLogLog.png"}],msdPlotLogLog];*)

msd=N[Mean[msdData[[All,2]]]];

{msd,msdExponent}
]


(* ::Subsection:: *)
(*calculateTemporalFreqSpectrum*)


(* ::Text:: *)
(*FFT with windowing*)
(*see: https://en.wikipedia.org/wiki/Window_function#Blackman_window*)
(**)
(*freqs f = [0,0.5*nSamples]/(dt*nSamples)*)


(* ::Input:: *)
(*{freqs,powerSpectralDensity,phases}=calculateTemporalFreqSpectrum[dt,signal,windowingQ,nPaddedSamples]*)


(* ::Subsubsection::Closed:: *)
(*test: frequency identification of n full periods*)


(* ::Input:: *)
(*dt=0.01;*)
(*f=2.0;*)
(*nSamples=1000;*)
(*times=Range[0,nSamples]dt;*)
(*sData2=Table[Sin[2\[Pi] f dt step](*+Sin[2\[Pi] 1.1f dt step]*),{step,0,nSamples,1}];*)
(*{windowingQ,nPadding}={True,Length[sData2]};*)
(*plotSignal=ListLinePlot[{times,sData2}\[Transpose],(*Evaluate[ops],*)Frame->True,FrameLabel->{"time t (s)","state s"},FrameStyle->Directive[Black,20,AbsoluteThickness[2]],ImageSize->400];*)
(*psdIdeal=calculateTemporalFreqSpectrum[dt,sData2,windowingQ,nPadding][[{1,2}]]\[Transpose](*\[LeftDoubleBracket];;100\[RightDoubleBracket]*);*)
(*plotPSD=ListLinePlot[psdIdeal,FrameLabel->{"frequency f (Hz)","\!\(\*SubscriptBox[\(PSD\), \(s\)]\)"},PlotRange->{{0,10},(*{0,0.4}*)All},Joined->True,ImageSize->400,Frame->True,FrameStyle->Directive[Black,20,AbsoluteThickness[2]],Mesh->All];*)
(*Grid[{{plotSignal,plotPSD}}]*)
(**)
(*{windowingQ,nPadding}={True,10Length[sData2]};*)
(*f=0.1;*)
(*times=Range[0,nSamples]dt;*)
(*sData2=Table[Sin[2\[Pi] f dt t]+Sin[2\[Pi] 1.1f dt t],{t,0,nSamples,1}];*)
(*plotSignal=ListLinePlot[{times,sData2}\[Transpose],(*Evaluate[ops],*)Frame->True,FrameLabel->{"time t (s)","state s"}]*)
(*psdIdeal=calculateTemporalFreqSpectrum[dt,sData2,windowingQ,nPadding][[{1,2}]]\[Transpose](*\[LeftDoubleBracket];;100\[RightDoubleBracket]*);*)
(*plotPSD=ListPlot[psdIdeal,FrameLabel->{"frequency f (Hz)","PSD(s)"},PlotRange->{{0,0.2},{0,0.4}(*All*)},Joined->True,ImageSize->400,Frame->True,FrameStyle->Directive[Black,20],Mesh->All]*)


(* ::Subsubsection::Closed:: *)
(*test: frequency identification of q partial periods*)


(* ::Input:: *)
(*f=0.1;*)
(*dt=0.01;*)
(*nSamples=9700;*)
(*times=Range[0,nSamples]dt;*)
(*sData2=Table[Sin[2\[Pi] f dt t],{t,0,nSamples,1}];*)
(*ps=ListLinePlot[{times,sData2}\[Transpose],(*Evaluate[ops],*)FrameLabel->{"time t (s)","state s"}]*)
(*psdNaive=calculateTemporalFreqSpectrum[dt,sData2,False][[{1,2}]]\[Transpose][[;;600]];*)
(*psdCorrected=calculateTemporalFreqSpectrum[dt,sData2,True][[{1,2}]]\[Transpose][[;;600]];*)
(*ListLogPlot[{psdIdeal,psdNaive,psdCorrected},FrameLabel->{"frequency f (Hz)","PSD(s)"},PlotRange->{{0,0.2},All},Joined->True,ImageSize->400,Frame->True,FrameStyle->Directive[Black,20]]*)


(* ::Subsubsection::Closed:: *)
(*ftFrequencyIdentification (deprecated)*)


(* ::Code:: *)
(*ftFrequencyIdentification[timeseries_,timeseriesDuration_]:=Module[{nSamples,dt,spectrumLength,powerSpectralDensity,freqs,peakThreshold,foundFrequency},nSamples=Length@timeseries;*)
(*dt=timeseriesDuration/nSamples;*)
(*Print["sampling frequency/rate = ",1.0/dt,"\nNyquist freq: ",(1.0/dt-1)*0.5];*)
(*spectrumLength=Floor[0.5 nSamples];*)
(*powerSpectralDensity=Abs[Fourier[timeseries][[;;spectrumLength]]]^2;*)
(*freqs=Range[0,spectrumLength-1]/(dt*nSamples);*)
(*peakThreshold=0.1;*)
(*foundFrequency=freqs[[#]]&/@Round[FindPeaks[powerSpectralDensity,0.0,0,peakThreshold][[All,1]]];*)
(*Print["found frequency: ",foundFrequency];*)
(*{foundFrequency,freqs,powerSpectralDensity}]*)


(* ::Subsubsection:: *)
(*analyze data*)


(* ::Input:: *)
(*data=Import[FileNameJoin[{$HomeDirectory,"Desktop","intensityTimetracesForJan.mat"}]][[1,1]];*)
(*(** Assuming all time points have the same time difference between them **)*)
(*dt=Mean[Differences@data[[1,All,1]]]*)
(*values=data[[All,All,2]];*)
(**)
(*(** fourier settings **)*)
(*windowingQ=False;*)
(*nPadding=2Length[values[[1]]];*)
(*padMemory=ConstantArray[0.0,nPadding];*)
(**)
(*(** calculate ffts **)*)
(*ffts=calculateTemporalFreqSpectrum[dt,#,windowingQ,nPadding,padMemory][[{1,2}]]&/@values;*)
(*freqData=ffts[[1,1]];*)
(*mPSDi=Mean[ffts[[All,2]]];*)
(*fMax=1.0;*)
(*fMaxIndex=FirstPosition[freqData,_?(#>=fMax&)][[1]];*)
(*mPSDiData={freqData,#/(Max[#]1.0)}\[Transpose][[2;;fMaxIndex]]&@mPSDi;*)
(**)
(*(** plot: intensity spectrum (PSD)**)*)
(*fTicks=Table[If[Mod[f,0.5]==0,{f,DecimalForm[f,{3,1}],{0.02,0}},{f,"",{0.01,0}}],{f,0,1,0.1}];*)
(*fTicks2={#[[1]],"",#[[3]]}&/@fTicks;*)
(*pTicks=Table[{t,DecimalForm[t,{2,1}]},{t,0,1,0.5}];*)
(*pTicks2={#[[1]],""}&/@pTicks;*)


(* ::Input:: *)
(*lineThickness=99.6923076923077`/72;*)
(*fontSize=10 99.6923076923077`/72;*)
(*plotOps={*)
(*		Axes->None*)
(*		,Frame->True*)
(*		,FrameStyle->Directive[Black,fontSize,AbsoluteThickness[lineThickness],FontFamily->"Arial"]*)
(*	};*)
(*ipx={54,12};*)
(*ListLinePlot[mPSDiData*)
(*	,Evaluate[plotOps]*)
(*	,AspectRatio->1/2*)
(*	,PlotRange->{{-0.01,1.01}fMax,{-0.05,1.25}}*)
(*	,FrameTicks->{{pTicks,pTicks2},{fTicks,fTicks2}}*)
(*	,ImagePadding->{ipx,{35,8}}*)
(*	,PlotStyle->Directive[Black,AbsoluteThickness[lineThickness]]*)
(*]*)


ListLinePlot[data[[1,All,{1,2}]]]


(* ::Subsubsection::Closed:: *)
(*debug*)


(* ::Input:: *)
(*Range[0,0.5 20000-1]/(dt 20000)//Max*)


(* ::Input:: *)
(*windowingQ=True;*)
(**)
(*(** original **)*)
(*signal=intensity;*)
(*nSamples=Length[signal];*)
(*psd1=Abs[Fourier[signal][[;;Floor[0.5nSamples]]]]^2;*)
(*psd1[[1]]*)
(*p1=ListLinePlot[signal,Frame->True];*)
(*p2=ListLogPlot[psd1,PlotRange->All,Joined->True,Frame->True];*)
(*Grid[{{p1,p2}}]*)
(**)
(*(** blackman original **)*)
(*signal1=If[windowingQ,*)
(*	a=0.16;*)
(*	windowFunctionBlackman=Table[(1-a)/2-1/2 Cos[2\[Pi] n/(nSamples-1)]+a/2 Cos[4\[Pi] n/(nSamples-1)],{n,0,nSamples-1}];*)
(*	signal windowFunctionBlackman,*)
(*	signal*)
(*];*)
(*ListPlot[windowFunctionBlackman];*)
(*psd2=Abs[Fourier[signal1][[;;Floor[0.5nSamples]]]]^2;*)
(*psd2[[1]]*)
(*p1=ListLinePlot[signal1,Frame->True];*)
(*p2=ListLogPlot[psd2,PlotRange->All,Joined->True,Frame->True];*)
(*Grid[{{p1,p2}}]*)
(**)
(*(** blackman zero-offset original **)*)
(*signal1b=signal-Mean[signal];*)
(*signal1c=If[windowingQ,*)
(*	a=0.16;*)
(*	windowFunctionBlackman=Table[(1-a)/2-1/2 Cos[2\[Pi] n/(nSamples-1)]+a/2 Cos[4\[Pi] n/(nSamples-1)],{n,0,nSamples-1}];*)
(*	signal1b windowFunctionBlackman,*)
(*	signal1b*)
(*];*)
(*psd2c=Abs[Fourier[signal1c][[;;Floor[0.5nSamples]]]]^2;*)
(*psd2c[[1]]*)
(*p1=ListLinePlot[signal1c,Frame->True];*)
(*p2=ListLogPlot[psd2c,PlotRange->All,Joined->True,Frame->True];*)
(*Grid[{{p1,p2}}]*)
(**)
(*(** zero-offset blackman original **)*)
(*signal2=signal1-Mean[signal1];*)
(*psd3=Abs[Fourier[signal2][[;;Floor[0.5nSamples]]]]^2;*)
(*psd3[[1]]*)
(*p1=ListLinePlot[signal2,Frame->True];*)
(*p2=ListLogPlot[psd3,PlotRange->All,Joined->True,Frame->True];*)
(*Grid[{{p1,p2}}]*)
(**)
(*(** zero-offset blackman zero-offset original **)*)
(*signal2b=signal1c-Mean[signal1c];*)
(*psd3b=Abs[Fourier[signal2b][[;;Floor[0.5nSamples]]]]^2;*)
(*psd3b[[1]]*)
(*p1=ListLinePlot[signal2b,Frame->True];*)
(*p2=ListLogPlot[psd3b,PlotRange->All,Joined->True,Frame->True];*)
(*Grid[{{p1,p2}}]*)


(* ::Input:: *)
(*ListLogPlot[Abs@(psd2-psd2c)[[2;;]],Joined->True,PlotRange->All];*)


(* ::Input:: *)
(*{iStart,iEnd}={2,100};*)
(*ListLogPlot[{psd1[[iStart;;iEnd]],psd2[[iStart;;iEnd]],psd2c[[iStart;;iEnd]],psd3[[iStart;;iEnd]],psd3b[[iStart;;iEnd]]},Joined->True,PlotLegends->{"f","window @ f","window @ mean @ f","mean @ window @ f","mean @ window @ mean @ f"},Frame->True,AspectRatio->1]*)


(* ::Subsubsection:: *)
(*code*)


calculateTemporalFreqSpectrum[dt_,signal_,windowingQ_:False,nPaddedSamples_:0,paddedMemory_:{}]:=Module[
{nSamples,signal0,signal1,signal2,signal3,a,windowFunctionBlackman,spectrumLength,powerSpectralDensity,freqs,ft,phases},
nSamples=Length[signal];

(** remove offset to avoid major peak at zeroth wavenumber **)
signal0=signal-Mean[signal];

(** window function since signal is usually non-periodic at boundaries **)
signal2=If[windowingQ,
	a=0.16;
	windowFunctionBlackman=Table[(1-a)/2-1/2 Cos[2\[Pi] n/(nSamples-1)]+a/2 Cos[4\[Pi] n/(nSamples-1)],{n,0,nSamples-1}];
	signal0 windowFunctionBlackman,
	signal0
];

(** padding for better resolution **)
If[nPaddedSamples>Length[signal],
	(*nPaddedSamples=3nSamples;*)
	signal3=paddedMemory;
	signal3[[;;Length[signal2]]]=signal2;  (** fast alternative to native padding function **)
	(*signal2=PadRight[signal1,nPaddedSamples];*) (** pad signal with zeros, so that it has total length of nPaddedSamples; slow! **)
	nSamples=Length[signal3];,
	signal3=signal2;
];

(** fourier transform **)
spectrumLength=Floor[0.5nSamples];
ft=Fourier[signal3,FourierParameters->{1,1}][[;;spectrumLength]];
phases=Arg[ft];
powerSpectralDensity=dt/nSamples Abs[ft]^2;
freqs=Range[0,spectrumLength-1]/(dt nSamples);

(** output **)
{freqs,powerSpectralDensity,phases}
]


(* ::Subsection::Closed:: *)
(*meanStd*)


meanStd[vals_]:={Mean@#,StandardDeviation@#}&@vals


(* ::Subsection::Closed:: *)
(*oneDriveSync*)


(* ::Text:: *)
(*start/stop one drive service*)
(** uses windows cmd, not powershell!*)


(* ::Input:: *)
(*oneDriveSync["stop"]*)


(* ::Input:: *)
(*oneDriveSync["start"]*)


(* ::Subsubsection::Closed:: *)
(*code*)


oneDriveSync[state_]:=Module[{},
	If[$OperatingSystem=="Windows",
		Switch[state
			,"start",
				Run["start \"OneDrive\" /B \"C:\\Program Files\\Microsoft OneDrive\\OneDrive.exe\" /background"];
			,"stop",
				Run["\"C:\\Program Files\\Microsoft OneDrive\\OneDrive.exe\" /shutdown"];
			,_,
				Print["Error: Unknown state ("<>state<>")!"];
		];
	,
		Print["Operating system is not supported!"];
	];
]


(* ::Subsection::Closed:: *)
(*phaseUnwrap*)


(* ::Text:: *)
(*input: list of real phases*)


phaseUnwrap=Compile[{{args,_Real,1}},Module[{},args+Prepend[2 Pi Accumulate@-IntegerPart@Differences[args/Pi],0]],RuntimeOptions->"EvaluateSymbolically"->False];


(* ::Subsection::Closed:: *)
(*progressIndicatorRemainingTime, updateProgressIndicator*)


(* ::Text:: *)
(*red warnings in dynamic context are erroneous*)


(* ::Subsubsection:: *)
(*use case*)


(* ::Input:: *)
(*progressIndicatorRemainingTime[];*)
(**)
(*nIterations=100;*)
(*Do[*)
(*Pause[0.02];*)
(*updateProgressIndicator[nIterations];*)
(*,{i,nIterations}];*)


(* ::Subsubsection:: *)
(*code*)


progressIndicatorRemainingTime[]:=Module[{(*startTime,*)endTime,remainingTime},

(** local variables **)
startTime=AbsoluteTime[];

(** global variables **)
counter=0;
fractionDone=0.0;
currentTime=startTime;

(** dynamic indicator **)
PrintTemporary[Overlay[{ProgressIndicator[Dynamic[fractionDone],{0,1}],Dynamic[
endTime=If[fractionDone>0,(currentTime-startTime)/fractionDone,0.0];
remainingTime=endTime-(currentTime-startTime);
ToString[DecimalForm[100.0 fractionDone,{4,1}]]<>"\[ThinSpace]% | "<>IntegerString[Round[remainingTime]]<>"\[ThinSpace]s left"]},Alignment->"Center"]];
]


updateProgressIndicator[nIterations_]:=Module[{},
fractionDone=N[++counter/nIterations];
currentTime=AbsoluteTime[];
]


(* ::Subsection::Closed:: *)
(*renderVideo*)


(* ::Text:: *)
(*with ffmpeg*)
(*TODO: throw error if ffmpeg is not installed*)


(* ::Input:: *)
(*renderVideo[pthInput,pthFnOutput,fps];*)


(* ::Subsubsection:: *)
(*debug*)


(* ::Input:: *)
(*pthInput=FileNameJoin[{$HomeDirectory,"Desktop","mmamovie_intensity_images_young"}];*)
(*pthFnOutput=FileNameJoin[{pthVideo,"Supplementary Videos","Supplementary Video 9 intensity images.mp4"}];*)
(*nImages=Length[FileNames[All,pthInput]];*)
(*fps=50;*)
(*crf=10;*)
(*code="ffmpeg -y -framerate "<>ToString[fps]<>" -i \""<>FileNameJoin[{pthInput,"p_%0"<>ToString[IntegerLength[nImages]+1]<>"d.png"}]<>"\" -vf scale=\"trunc(iw/2)*2:trunc(ih/2)*2\" -vcodec libx264 -pix_fmt yuv420p -preset slow -crf "<>ToString[crf]<>" -r "<>ToString[fps]<>" \""<>pthFnOutput<>"\""*)


(* ::Input:: *)
(*Run[code]*)


(* ::Subsubsection:: *)
(*code*)


renderVideo[pthInput_,pthFnOutput_,fps_]:=Module[{imageFiles,nDigits,crf,code},
	imageFiles=FileNames[All,pthInput];
	nDigits=ToString[StringLength[StringSplit[FileBaseName[imageFiles[[1]]],"_"][[2]]]];
	crf=10;
	code="ffmpeg -y -framerate "<>ToString[fps]<>" -i \""<>FileNameJoin[{pthInput,"p_%0"<>nDigits<>"d.png"}]<>"\" -vf scale=\"trunc(iw/2)*2:trunc(ih/2)*2\" -vcodec libx264 -pix_fmt yuv420p -preset slow -crf "<>ToString[crf]<>" -r "<>ToString[fps]<>" \""<>pthFnOutput<>"\"";
	Run[code];
	
	(** debug output **)
	code
]


(* ::Subsection::Closed:: *)
(*sphericalNorm*)


sphericalNorm[v_,\[Theta]_]:=v . DiagonalMatrix[{1,Sin[\[Theta]]^2}] . v


(* ::Subsection::Closed:: *)
(*spherical base and tangent vectors er, e\[Theta], e\[Phi], r\[Theta], r\[Phi]*)


er[\[Theta]_,\[Phi]_]:={Sin[\[Theta]]Cos[\[Phi]],Sin[\[Theta]]Sin[\[Phi]],Cos[\[Theta]]};
e\[Theta][\[Theta]_,\[Phi]_]:={Cos[\[Theta]]Cos[\[Phi]],Cos[\[Theta]]Sin[\[Phi]],-Sin[\[Theta]]};
e\[Phi][\[Theta]_,\[Phi]_]:={-Sin[\[Phi]],Cos[\[Phi]],0};
r\[Theta][\[Theta]_,\[Phi]_]:={Cos[\[Theta]]Cos[\[Phi]],Cos[\[Theta]]Sin[\[Phi]],-Sin[\[Theta]]};
r\[Phi][\[Theta]_,\[Phi]_]:={-Sin[\[Theta]]Sin[\[Phi]],Sin[\[Theta]]Cos[\[Phi]],0};


(* ::Section::Closed:: *)
(*special analysis*)


(* ::Subsection::Closed:: *)
(*avgFreqThetaNeuron*)


avgFreqThetaNeuron[u0Dim_]:=Sqrt[(((rDropletDim+rEcoliDim)^2 u0Dim^2 dragSwimDim^2-gDim^2 mDim^2 rcomDim^2)/((rDropletDim+rEcoliDim)^4 dragDim^2))]/(2.0\[Pi])


(* ::Subsection::Closed:: *)
(*fGrav*)


fGrav=m gDim rcomDim/(rDropletDim+rEcoliShortDim){0,0,1};


(* ::Subsection::Closed:: *)
(*fSwim*)


fSwim[p_]:=swimFactor dragTransDim u0Dim p;


(* ::Subsection::Closed:: *)
(*freqOnlyRunning*)


freqOnlyRunning[u0Dim_]:=swimFactor dragTransDim u0Dim/(2\[Pi] (rDropletDim+rEcoliDim));


(* ::Subsection::Closed:: *)
(*maxTheta*)


Clear[maxTheta];
maxTheta[u0Dim_,fudgeFactor_:1.0]:=ArcSin[((rDropletDim+rEcoliShortDim) dragSwimDim u0Dim)/(fudgeFactor mDim gDim rcomDim)]


(* ::Section::Closed:: *)
(*droplet geometry*)


(* ::Subsection::Closed:: *)
(*calcd*)


(* ::Text:: *)
(*Note: For vr!=1, ri can be < rd -> select criterion fails atm!*)


(* ::Input:: *)
(*{rd,ri,vr}={1.0,100.0,1.0};*)
(*d=calcd[rd,ri,vr]*)


(* ::Subsubsection::Closed:: *)
(*test*)


(* ::Input:: *)
(*{rd,ri,vr}={1.0,100.0,1.0};*)
(*d=calcd[rd,ri,vr]*)


(* ::Input:: *)
(*{rd,ri,vr}={3.5*^-6,100.0*^-6,1.0};*)
(*d=calcd[rd,ri,vr]*)


(* ::Subsubsection::Closed:: *)
(*C++ debug*)


(* ::Input:: *)
(*Clear[rd,ri,vr]*)
(*(dd+rd-ri)^2  (dd^2 + 2dd(ri-rd) - 3(rd+ri)^2) + 16dd rd^3/(1+vr)==0//Expand(*Expand[#,dd]&*)*)
(*rd=1.0;*)
(*ri=4.0;*)
(*vr=1.0;*)
(*(dd+rd-ri)^2  (dd^2 + 2dd(ri-rd) - 3(rd+ri)^2) + 16dd rd^3/(1+vr)==0//Expand*)
(*{-3 rd^4-3 ri^4+6 rd^2 ri^2,-8 rd^3+8 ri^3+(16 rd^3)/(1+vr),-6(rd^2+ ri^2),0.0,1.0}*)
(*NSolve[(dd+rd-ri)^2  (dd^2 + 2dd(ri-rd) - 3(rd+ri)^2) + 16dd rd^3/(1+vr)==0]*)


(* ::Subsubsection:: *)
(*code*)


calcd[rd_,ri_,vr_]:=Module[{ds0},
Clear[dd];
ds0=dd/.NSolve[(dd+rd-ri)^2  (dd^2 + 2dd(ri-rd) - 3(rd+ri)^2) + 16dd rd^3/(1+vr)==0,dd,Reals];
If[#!={},#[[-1]],Nothing]&@Select[ds0,(-rd<#-ri<rd)&]
];


(* ::Subsection::Closed:: *)
(*calczcom*)


(* ::Text:: *)
(*source: Sara thesis eq A.30*)


(* ::Subsubsection::Closed:: *)
(*test*)


(* ::Input:: *)
(*Module[{\[Rho]FDim=1.01*^3,\[Rho]HDim=1.44*^3,volumeRatio=1.0,rDropletDim=3.5*^-6,rInterface=100.0*^-6,d,m,zcom},*)
(*d=calcd[rDropletDim,rInterface,volumeRatio];*)
(*m=mass[volumeRatio,\[Rho]FDim,\[Rho]HDim,rDropletDim];*)
(*zcom=calczcom[rDropletDim,rInterface,d,\[Rho]FDim,\[Rho]HDim,m]*)
(*]*)


(* ::Input:: *)
(*Module[{\[Rho]FDim=1.01*^3,\[Rho]HDim=1.44*^3,volumeRatio=1.0,rDropletDim=3.5*^-6,rInterface=3.5*^-6,d,m,zcom},*)
(*d=calcd[rDropletDim,rInterface,volumeRatio];*)
(*m=mass[volumeRatio,\[Rho]FDim,\[Rho]HDim,rDropletDim];*)
(*zcom=calczcom[rDropletDim,rInterface,d,\[Rho]FDim,\[Rho]HDim,m]*)
(*]*)


(* ::Subsubsection::Closed:: *)
(*code*)


ii[rd_,ri_,d_]:=(rd^2-ri^2+d^2)/(2.0d);
calczcom[rd_,ri_,d_,\[Rho]F_,\[Rho]H_,m_]:=(\[Rho]H-\[Rho]F)/(12m)\[Pi]( (d-ri)^2 * (d^2-3ri^2+2d ri) + 3rd^4 - 4 d ii[rd,ri,d]^3 );


(* ::Subsection::Closed:: *)
(*mass*)


(* ::Input:: *)
(*m=mass[volumeRatio,\[Rho]FDim,\[Rho]HDim,rDropletDim]*)


(* ::Subsubsection:: *)
(*test*)


(* ::Input:: *)
(*Module[{\[Rho]FDim=1.01*^3,\[Rho]HDim=1.44*^3,volumeRatio=1.0,rDropletDim=3.5*^-6},mass[volumeRatio,\[Rho]FDim,\[Rho]HDim,rDropletDim]]*)


(* ::Input:: *)
(*2.200031155370152`*^-13*10^13*)


(* ::Subsubsection:: *)
(*code*)


(** A.17 Sara's thesis **)
mass[vr_,\[Rho]F_,\[Rho]H_,rd_]:=4.0/3.0*(\[Rho]F+vr \[Rho]H)/(1.0+vr)\[Pi] rd^3;


(* ::Section::Closed:: *)
(*projectors*)


(* ::Subsection::Closed:: *)
(*pCircle*)


pCircle[r_,\[Theta]_]:=Simplify[sphereSurfaceProjector2d[r{Sin[\[Theta]],Cos[\[Theta]]},{0,0}],r>0\[And]\[Theta]\[Element]Reals]


(* ::Subsection::Closed:: *)
(*pSphere*)


pSphere[r_,\[Theta]_,\[Phi]_]:=Simplify[sphereSurfaceProjector[r{Cos[\[Phi]]Sin[\[Theta]],Sin[\[Phi]]Sin[\[Theta]],Cos[\[Theta]]},{0,0,0}],r>0\[And]\[Theta]\[Element]Reals\[And]\[Phi]\[Element]Reals]


(* ::Subsection::Closed:: *)
(*sphereSurfaceProjector2d*)


sphereSurfaceProjector2d[{x_,y_},circleOrigin_]:=Module[{n},
n=Normalize[{x,y}-circleOrigin];
IdentityMatrix[Length@n]-KroneckerProduct[n,n]
]


(* ::Subsection::Closed:: *)
(*sphereSurfaceProjector*)


(* ::Subsubsection::Closed:: *)
(*test*)


(* ::Input:: *)
(*sphereSurfaceProjector[er[\[Theta],\[Phi]]]//ComplexExpand//Simplify*)


(* ::Subsubsection::Closed:: *)
(*code*)


sphereSurfaceProjector[{x_,y_,z_}]:=Module[{n},
	n=Normalize[{x,y,z}];
	IdentityMatrix[Length@n]-KroneckerProduct[n,n]
]

sphereSurfaceProjector[{x_,y_,z_},sphereOrigin_]:=Module[{n},
n=Normalize[{x,y,z}-sphereOrigin];
IdentityMatrix[Length@n]-KroneckerProduct[n,n]
]


(* ::Subsection::Closed:: *)
(*surfaceProjection*)


surfaceProjection[v_,nUnit_]:=(1-TensorProduct[#,#]&@nUnit) . v


(* ::Section::Closed:: *)
(*utilities*)


(* ::Subsection::Closed:: *)
(*circle3D*)


circle3D[centre_:{0,0,0}, radius_:1, normal_:{0,0,1}, angle_:{0, 2 Pi}]:=
  Composition[
    Line,
    Map[RotationTransform[{{0, 0, 1}, normal}, centre], #] &,
    Map[Append[#, Last@centre] &, #] &,
    Append[DeleteDuplicates[Most@#], Last@#] &,
    Level[#, {-2}] &,
    MeshPrimitives[#, 1] &,
    DiscretizeRegion,
    If
  ][
    First@Differences@angle >= 2 Pi,
    Circle[Most@centre, radius],
    Circle[Most@centre, radius, angle]
  ]


(* ::Subsection::Closed:: *)
(*colorbar*)


(* ::Input:: *)
(*Show[colorbar[{0,1},grayTicks,False,fs,is,"test",2,GrayLevel],ImagePadding->{{1,35},{60,10}}]*)


colorbar[{min_,max_},ticks_,leftQ_,fs_,is_,label_,absThick_,colorFunction_: Automatic, rotateLabelQ_:True, divs_: 150]:=
DensityPlot[y,{x,0,0.1},{y,min,max},
AspectRatio->30,
ColorFunction->colorFunction,
ColorFunctionScaling->(*False*)True,
Frame->True,
FrameLabel->{{None,label},{None,None}},
FrameStyle->Directive[Black,AbsoluteThickness[absThick]],
FrameTicks->{If[leftQ,{{#[[1]],#[[2]],{0,0.5}}&/@ticks,None},{None,{#[[1]],#[[2]],{0,0.5}}&/@ticks}],{None,None}},
FrameTicksStyle->Directive[FontFamily->"Helvetica",fs,Black],
ImageSize->{Automatic,is},
LabelStyle->Directive[FontFamily->"Helvetica",fs,Black],
MaxRecursion->0,
PlotRangePadding->0,
PlotPoints->{2,divs},
RotateLabel->rotateLabelQ
]


(* ::Subsection::Closed:: *)
(*log10Subdivide*)


log10Subdivide[start_?Positive,end_?Positive,increments_]:=10^Range[Log10@start,Log10@end,Log10[end/start]/increments]


(* ::Subsection::Closed:: *)
(*progressIndicatorRemainingTime, updateProgressIndicator*)


(* ::Text:: *)
(** for interactive execution of code*)
(** red warnings in dynamic context are erroneous*)


(* ::Subsubsection::Closed:: *)
(*use case*)


(* ::Input:: *)
(*progressIndicatorRemainingTime[];*)
(**)
(*nIterations=100;*)
(*Do[*)
(*Pause[0.02];*)
(*updateProgressIndicator[nIterations];*)
(*,{i,nIterations}];*)


(* ::Subsubsection:: *)
(*code*)


progressIndicatorRemainingTime[]:=Module[{startTime,endTime,remainingTime},

(** local variables **)
startTime=AbsoluteTime[];

(** global variables **)
counter=0;
fractionDone=0.0;
currentTime=startTime;

(** dynamic indicator **)
PrintTemporary[Overlay[{ProgressIndicator[Dynamic[fractionDone],{0,1}],Dynamic[
endTime=If[fractionDone>0,(currentTime-startTime)/fractionDone,0.0];
remainingTime=endTime-(currentTime-startTime);
ToString[DecimalForm[100.0 fractionDone,{4,1}]]<>"\[ThinSpace]% | "<>IntegerString[Round[remainingTime]]<>"\[ThinSpace]s left"]},Alignment->"Center"]];
]


updateProgressIndicator[nIterations_]:=Module[{},
fractionDone=N[++counter/nIterations];
currentTime=AbsoluteTime[];
]


(* ::Chapter:: *)
(*Simulation scan functions*)


(* ::Subsection::Closed:: *)
(*analyzeHeterogeneityData*)


analyzeHeterogeneityData[u0Dim_,hetQ_,pthout_]:=Module[{pthFileName,fstream,type,heterogeneityData,mean,std,plotLabel,plot},
If[hetQ!=0,
	pthFileName=FileNameJoin[{NotebookDirectory[],"build","u0_heterogeneity.bin"}];
	fstream=OpenRead[pthFileName,BinaryFormat->True];
	type="Real64";
	heterogeneityData=BinaryReadList[fstream,type];  (** read data **)
	Close[fstream];
	
	(** rescale data **)
	heterogeneityData*=10^6 dragDim rDropletDim/dragSwimDim;    (** in \[Mu]m/s instead of m/s! invert C++ program scaling: activeFrictionT*u0/(rDropletDim*dragFactor) **)
	{mean,std}=ToString[DecimalForm[#,{3,2}]]&/@{Mean@#,StandardDeviation@#}&@heterogeneityData;
	plotLabel="\!\(\*SubscriptBox[\(u\), \(0\)]\) ="<>ToString[DecimalForm[10^6 u0Dim,{3,1}]]<>" = "<>mean<>" \[PlusMinus] "<>std<>" \[Mu]m/s";
	plot=Show[
		Histogram[heterogeneityData,{0,30,0.5},Frame->True,FrameLabel->{"speed \!\(\*
StyleBox[SubscriptBox[\"u\", 
RowBox[{\"0\", \" \"}]],\nFontSlant->\"Italic\"]\)(\[Mu]m/s)","occurrence"},
		FrameStyle->Directive[Black,18,AbsoluteThickness[2]],AspectRatio->1,PlotLabel->Style[plotLabel,Black,20],
		Epilog->{Black,InfiniteLine[{{#,0},{#,1}}]&@(u0Dim 10^6)}]
	];
	Export[FileNameJoin[{pthout,"u0_heterogeneity_"<>ToString[DecimalForm[10^6 u0Dim,{3,1}]]<>".png"}],plot];
	Print[plot];
];
]


(* ::Subsection::Closed:: *)
(*calcFits*)


(* ::Subsubsection:: *)
(*debug*)


(* ::Input:: *)
(*fMax=1.0;*)
(*fMaxIndex=FirstPosition[freqs,_?(#>=fMax&)][[1]];*)
(*{indexStart,indexEnd}={2,fMaxIndex};*)
(*iFitTheory=NonlinearModelFit[{freqs,#/Max[#]}\[Transpose][[indexStart;;indexEnd]],{aa/(1+(\[Omega]-\[Omega]0)^2/bb^2),aa>0,bb>0,\[Omega]0>0},{{aa,1.0},{bb,0.01},{\[Omega]0,0.001}},\[Omega]]&@mPSDi;*)


(* ::Input:: *)
(*iFitTheory=NonlinearModelFit[{freqs,mPSD\[Theta]}\[Transpose][[indexStart;;indexEnd]],{aa/(1+(\[Omega]-\[Omega]0)^2/bb^2),aa>0,bb>0,\[Omega]0>0},{{aa,1.0},{bb,0.01},{\[Omega]0,0.001}},\[Omega]];*)


(* ::Subsubsection:: *)
(*code*)


calcFits[fMax_,freqs_,mPSDs_,mPSD\[Theta]_,mPSDi_]:=Module[{complexFitQ,fMaxIndex,indexStart,indexEnd,sFit,\[Theta]Fit,iFit,\[Theta]FitTheory,weights,iFitTheory,sWidth,\[Theta]Width,iWidth,meanThetaFit},

Clear[\[Omega],\[Omega]0,aa,bb,cc,i0];
fMaxIndex=FirstPosition[freqs,_?(#>=fMax&)][[1]];
{indexStart,indexEnd}={2,fMaxIndex};

(** fitting normalized average spectrum in s,\[Theta],I **)
(*,Method\[Rule]{NMinimize}*)
(*{sFit,\[Theta]Fit,iFit}=NonlinearModelFit[{freqs,#/Max[#]}\[Transpose]\[LeftDoubleBracket]indexStart;;indexEnd\[RightDoubleBracket],{aa/(1+(\[Omega]-\[Omega]0)^2/bb^2),aa>0,bb>0,\[Omega]0>=0},{{aa,1},{bb,0.1},{\[Omega]0,0.001}},\[Omega],MaxIterations->1000]&/@{mPSDs,mPSD\[Theta],mPSDi};*)
(** Lorentz with peak at 0 Hz **)
{sFit,\[Theta]Fit,iFit}=NonlinearModelFit[{freqs,#/Max[#]}\[Transpose][[indexStart;;indexEnd]],{aa/(1+(\[Omega]/bb)^2),aa>0,bb>0},{{aa,1},{bb,0.1}},\[Omega],MaxIterations->1000]&/@{mPSDs,mPSD\[Theta],mPSDi};
{sWidth,\[Theta]Width,iWidth}=Abs[#["BestFitParameters"][[2,2]]]&/@{sFit,\[Theta]Fit,iFit};
{\[Theta]FitTheory,iFitTheory}={\[Theta]Fit,iFit};
meanThetaFit=0.0;

(** TODO: better control structure; some kind of switch **)
complexFitQ=False;
If[complexFitQ,
\[Theta]FitTheory=NonlinearModelFit[{freqs,#/Max[#]}\[Transpose][[indexStart;;indexEnd]],{aa/(1+(\[Omega]-\[Omega]0)^2/bb^2),aa>0,bb>0,\[Omega]0>=0},{{aa,1.0},{bb,0.1},{\[Omega]0,0.001}},\[Omega]]&@mPSD\[Theta];
iFitTheory=NonlinearModelFit[{freqs,#/Max[#]}\[Transpose][[indexStart;;indexEnd]],{aa/(1+(\[Omega]-\[Omega]0)^2/bb^2),aa>0,bb>0,\[Omega]0>=0},{{aa,1.0},{bb,0.01},{\[Omega]0,0.001}},\[Omega]]&@mPSDi;
(*iFitTheory=NonlinearModelFit[{freqs,#}\[Transpose]\[LeftDoubleBracket]indexStart;;indexEnd\[RightDoubleBracket],{i0^2 (\[Pi]^2/4+(2 aa^2 bb^3 Sqrt[2/\[Pi]])/(4 bb^2+(\[Omega]-2 \[Omega]0)^2)+(aa bb^2(2  cc Sqrt[2/\[Pi]]- \[Pi]))/(bb^2+(\[Omega]-\[Omega]0)^2))^2,aa>0,bb>0,cc>0,\[Omega]0>0},{{i0,0.01},{aa,\[Theta]FitTheory["BestFitParameters"]\[LeftDoubleBracket]1,2\[RightDoubleBracket]},{bb,\[Theta]FitTheory["BestFitParameters"]\[LeftDoubleBracket]2,2\[RightDoubleBracket]},{cc,1.01},{\[Omega]0,\[Theta]FitTheory["BestFitParameters"]\[LeftDoubleBracket]3,2\[RightDoubleBracket]}},\[Omega]]&@mPSDi;*)
weights=ConstantArray[1.0,indexEnd-indexStart+1];
(*iFitTheory=NonlinearModelFit[{freqs,mPSDi}\[Transpose]\[LeftDoubleBracket]indexStart;;indexEnd\[RightDoubleBracket],{aa bb^2i0^2 ((2 aa bb Sqrt[2/\[Pi]])/(4 bb^2+(\[Omega]-2 \[Omega]0)^2)+(2  cc Sqrt[2/\[Pi]]- \[Pi])/(bb^2+(\[Omega]-\[Omega]0)^2))^2/.\[Theta]FitTheory["BestFitParameters"],cc>0},{{i0,0.1},{cc,2.1}},\[Omega],Weights\[Rule]weights];*)
(*iFitTheory=NonlinearModelFit[{freqs,mPSDi}\[Transpose]\[LeftDoubleBracket]indexStart;;indexEnd\[RightDoubleBracket],{i0^2 (1/(2520 Sqrt[2 \[Pi]])aa bb^2 (-((128 aa^7 bb^7)/(64 bb^2+(\[Omega]-8 \[Omega]0)^2))+(448 aa^6 bb^6 (-2 cc+\[Pi]))/(49 bb^2+(\[Omega]-7 \[Omega]0)^2)-(672 aa^5 bb^5 (-2+4 cc^2-4 cc \[Pi]+\[Pi]^2))/(36 bb^2+(\[Omega]-6 \[Omega]0)^2)+(560 aa^4 bb^4 (-8 cc^3+12 cc^2 \[Pi]+\[Pi] (-6+\[Pi]^2)-6 cc (-2+\[Pi]^2)))/(25 bb^2+(\[Omega]-5 \[Omega]0)^2)-(280 aa^3 bb^3 (24+16 cc^4-32 cc^3 \[Pi]-12 \[Pi]^2+\[Pi]^4-8 cc \[Pi] (-6+\[Pi]^2)+24 cc^2 (-2+\[Pi]^2)))/(16 bb^2+(\[Omega]-4 \[Omega]0)^2)-(84 aa^2 bb^2 (32 cc^5-80 cc^4 \[Pi]-40 cc^2 \[Pi] (-6+\[Pi]^2)+80 cc^3 (-2+\[Pi]^2)-\[Pi] (120-20 \[Pi]^2+\[Pi]^4)+10 cc (24-12 \[Pi]^2+\[Pi]^4)))/(9 bb^2+(\[Omega]-3 \[Omega]0)^2)-(14 aa bb (-720+64 cc^6-192 cc^5 \[Pi]+360 \[Pi]^2-30 \[Pi]^4+\[Pi]^6-160 cc^3 \[Pi] (-6+\[Pi]^2)+240 cc^4 (-2+\[Pi]^2)-12 cc \[Pi] (120-20 \[Pi]^2+\[Pi]^4)+60 cc^2 (24-12 \[Pi]^2+\[Pi]^4)))/(4 bb^2+(\[Omega]-2 \[Omega]0)^2)+1/(bb^2+(\[Omega]-\[Omega]0)^2)(-128 cc^7+448 cc^6 \[Pi]+560 cc^4 \[Pi] (-6+\[Pi]^2)-672 cc^5 (-2+\[Pi]^2)+84 cc^2 \[Pi] (120-20 \[Pi]^2+\[Pi]^4)-280 cc^3 (24-12 \[Pi]^2+\[Pi]^4)+\[Pi] (-5040+840 \[Pi]^2-42 \[Pi]^4+\[Pi]^6)-14 cc (-720+360 \[Pi]^2-30 \[Pi]^4+\[Pi]^6))))^2/.\[Theta]FitTheory["BestFitParameters"],\[Pi]/4<cc<\[Pi]},{{i0,0.1},{cc,2.1}},\[Omega],Weights\[Rule]weights];*)
(*iFitTheory=NonlinearModelFit[{freqs,mPSDi}\[Transpose]\[LeftDoubleBracket]indexStart;;indexEnd\[RightDoubleBracket],{i0^2 (\[Pi]^2/4+(2 aa^2 bb^3 Sqrt[2/\[Pi]])/(4 bb^2+(\[Omega]-2 \[Omega]0)^2)+(aa bb^2(2  cc Sqrt[2/\[Pi]]- \[Pi]))/(bb^2+(\[Omega]-\[Omega]0)^2))^2/.\[Theta]FitTheory["BestFitParameters"]\[LeftDoubleBracket]3\[RightDoubleBracket],cc>0},{{i0,0.0001},{aa,100\[Theta]FitTheory["BestFitParameters"]\[LeftDoubleBracket]1,2\[RightDoubleBracket]},{bb,\[Theta]FitTheory["BestFitParameters"]\[LeftDoubleBracket]2,2\[RightDoubleBracket]},{cc,10.01}},\[Omega]];*)
Print[\[Theta]FitTheory["BestFitParameters"]];
meanThetaFit=cc/.iFitTheory["BestFitParameters"];
];
Print[iFitTheory["BestFitParameters"]];


{sFit,\[Theta]Fit,iFit,\[Theta]FitTheory,iFitTheory,sWidth,\[Theta]Width,iWidth,meanThetaFit}
]


(* ::Subsection::Closed:: *)
(*calcForceFactors*)


(* ::Text:: *)
(*Calculate factors for spatially rescaled system *)
(*-> later: multiply space by rDropletDim to get out right values*)


(* ::Subsubsection::Closed:: *)
(*code*)


calcForceFactors[rDropletDim_,u0Dim_,swimMultiplier_:1.0]:=Module[{muDim,rEcoliDim,dragFactorT,gDim,vr,\[Rho]FDim,\[Rho]HDim,rEcoli,g,rDroplet,u0,mu,\[Rho]F,\[Rho]H,m,d,rcom,(*mDim,*)dDim,(*rcomDim,*)rRatio,swimFactor,gravFactor},
  muDim = 8.9*^-4;
  rEcoliDim = 0.5*^-6;
  gDim=9.81;
  vr=1.0;
  {\[Rho]FDim,\[Rho]HDim}={1.01*^3,1.44*^3};
  {rEcoli,g,rDroplet,u0}={rEcoliDim,gDim,rDropletDim,u0Dim}/rDropletDim;
  mu=muDim rDropletDim;
  {\[Rho]F,\[Rho]H}={\[Rho]FDim,\[Rho]HDim} rDropletDim^3;
  dragFactorT= 6.0 \[Pi] mu rEcoli;
  
  m=mass[vr,\[Rho]F,\[Rho]H,rDroplet];
  d=calcd[rDroplet,20rDroplet,vr];
  rcom=calczcom[rDroplet,20rDroplet,d,\[Rho]F,\[Rho]H,m];
  rRatio=rcom/(rDroplet+rEcoli);
  
  (** TODO: move through correct channels **)
  mDim=m;
  rcomDim=rcom rDropletDim;
  rEcoliShort=rEcoli;
  
  {swimFactor,gravFactor}={swimMultiplier dragFactorT u0, m g rRatio}
]


(* ::Subsection::Closed:: *)
(*calcGeometricalParameters*)


(* ::Input:: *)
(*{{rDropletDim,rEcoliShortDim,rEcoliLongDim,rInterfaceDim,\[Rho]FDim,\[Rho]HDim,dDim,mDim,rcomDim,volumeRatio},{rDroplet,rEcoliShort,rEcoliLong,rInterface,\[Rho]F,\[Rho]H,d,m,rcm}}=calcGeometricalParameters[]*)


calcGeometricalParameters[rDropletDim_:3.5*^-6]:=Module[{
rEcoliShortDim,rEcoliLongDim,rInterfaceDim,\[Rho]FDim,\[Rho]HDim,dDim,mDim,rcomDim,volumeRatio,
rDroplet,rEcoliShort,rEcoliLong,rInterface,\[Rho]F,\[Rho]H,d,m,rcm},
	
	rEcoliShortDim=0.5*^-6;
	rEcoliLongDim=1.0*^-6;
	rInterfaceDim=100.0*^-6;
	
	volumeRatio=1.0;
	
	{\[Rho]FDim,\[Rho]HDim}={1.01*^3,1.44*^3};
	
	(** in spatially rescaled units **)
	{rDroplet,rInterface,rEcoliShort,rEcoliLong}={rDropletDim,rInterfaceDim,rEcoliShortDim,rEcoliLongDim}/rDropletDim;
	{\[Rho]F,\[Rho]H}={\[Rho]FDim,\[Rho]HDim}rDropletDim^3;
	
	(** calculate derived quantities **)
	d=calcd[rDroplet,rInterface,volumeRatio];
	m=mass[volumeRatio,\[Rho]F,\[Rho]H,rDroplet];
	rcm=calczcom[rDroplet,rInterface,d,\[Rho]F,\[Rho]H,m];
	
	(** physical values **)
	mDim=m; (** Spatial units are "integrated" out. Only spatial units were rescaled. **)
	dDim=d rDropletDim;
	rcomDim=rcm  rDropletDim;
	
	(** output **)
	{{rDropletDim,rEcoliShortDim,rEcoliLongDim,rInterfaceDim,\[Rho]FDim,\[Rho]HDim,dDim,mDim,rcomDim,volumeRatio},{rDroplet,rEcoliShort,rEcoliLong,rInterface,\[Rho]F,\[Rho]H,d,m,rcm}}
]


(* ::Subsection::Closed:: *)
(*calcMeanPowerSpectra*)


(* ::Input:: *)
(*{freqs,mPSDs,mPSD\[Theta],mPSDi,stdPSD\[Theta],stdPSDi,sData,thetaData,phiData,iData,phasesS,phases\[Theta],phasesI,*)
(*meanThetas,meanMeanTheta,thetaNoise,mThetaNoiseMean,mThetaNoiseC,meanIminmax,histoTheta}=calcMeanPowerSpectra[trajectories,dt,gravFactor,plotNoiseQ,stepsTransient];*)


(* ::Text:: *)
(*Note: cutting off the initial transient makes the spectrum wider in case of otherwise flat trajectories*)


(* ::Subsubsection::Closed:: *)
(*code*)


calcMeanPowerSpectra[trajectories_,dt_,gravFactor_,plotNoiseQ_,stepsTransient_]:=Module[{windowingQ,(*psdData,*)sData,thetaData,phiData,iData,
psds,psd\[Theta],psdi,freqs,mPSDs,mPSD\[Theta],mPSDi,stdPSD\[Theta],stdPSDi,psd,phaseS,phase\[Theta],phaseI,phasesS,phases\[Theta],phasesI,meanTheta,meanMeanTheta,meanThetas,
thetaNoise,meanThetaNoise,mThetaNoiseMean,thetaNoiseC,mThetaNoiseC,nPadding,padMemory,thetaDataWrapped,phiDataWrapped,meanIminmax,histoTheta},

windowingQ=False;

(** get info from all simulation trajectories **)
nPadding=2Length[trajectories[[1]]];
padMemory=ConstantArray[0.0,nPadding];
histoTheta=ConstantArray[0.0,100];
psdData=Table[
	sData=trajectory[[All,7]];
	thetaData=ArcCos/@trajectory[[All,3]];   (** 0 - top pole (z=+1), \[Pi] - bottom pole (z=-1) **)
	histoTheta+=BinCounts[thetaData,{0,\[Pi],0.01\[Pi]}];
	If[plotNoiseQ,
		thetaNoise=thetaNoiseReconstruction[thetaData,dt,gravFactor];
		meanThetaNoise=Mean[thetaNoise];
		thetaNoiseC=CorrelationFunction[thetaNoise,{Length[thetaNoise]-1}];
	];
	phiData=phaseUnwrap[ArcTan[#[[1]],#[[2]]]&/@trajectory[[All,{1,2}]]];
	(** interpolation mirrored around \[Pi]/2 and \[Pi]/4 for \[Theta] and \[Phi] **)
	thetaDataWrapped=-Abs[thetaData-\[Pi]/2.0]+\[Pi]/2.0;
	phiDataWrapped=-Abs[Mod[phiData,\[Pi]/2.0]-\[Pi]/4.0]+\[Pi]/4.0;
	(*iData=Cos[#]^2 &/@ thetaData;*)
	(*iData=\[Theta]ii[#] &/@ thetaData;*)
	iData=MapThread[thetaPhiIntensityInterpolation[#1,#2]&,{thetaDataWrapped,phiDataWrapped}];
	(*iData=MapThread[fit\[Theta]\[Phi][#1,#2]&,{thetaData,phiData}];  *)
	meanTheta=Mean[thetaData];
	iMinMax=(Max[#]-Min[#])&@iData;
	{freqs,psd,phaseS}=calculateTemporalFreqSpectrum[dt,sData,windowingQ,nPadding,padMemory];
	psds={freqs,psd}\[Transpose];
	{freqs,psd,phase\[Theta]}=calculateTemporalFreqSpectrum[dt,thetaData,windowingQ,nPadding,padMemory];
	psd\[Theta]={freqs,psd}\[Transpose];
	{freqs,psd,phaseI}=calculateTemporalFreqSpectrum[dt,iData,windowingQ,nPadding,padMemory];
	psdi={freqs,psd}\[Transpose];
	{psds,psd\[Theta],psdi,phaseS,phase\[Theta],phaseI,meanThetaNoise,thetaNoiseC,meanTheta,iMinMax}
,{trajectory,trajectories[[All,stepsTransient;;]]}];
histoTheta/=100;

(** average info **)
mPSDs=Mean[psdData[[All,1,All,2]]];                      (** averaged power spectral density of noise s **)
mPSD\[Theta]=Mean[psdData[[All,2,All,2]]];                      (** averaged power spectral density of polar angle \[Theta] **)
mPSDi=Mean[psdData[[All,3,All,2]]];                      (** averaged power spectral density of intensity i **)
stdPSD\[Theta]=StandardDeviation[psdData[[All,2,All,2]]];       (** standard deviation of power spectral density of polar angle \[Theta] **)
stdPSDi=StandardDeviation[psdData[[All,3,All,2]]];       (** standard deviation of power spectral density of intensity i **)
phasesS=psdData[[All,4]];
phases\[Theta]=psdData[[All,5]];
phasesI=psdData[[All,6]];
{mThetaNoiseMean,mThetaNoiseC}=If[plotNoiseQ,{Mean[psdData[[All,7]]],Mean[psdData[[All,8]]]},{0,0}];
meanThetas=psdData[[All,9]];
meanMeanTheta=Mean[meanThetas];  (** scalar mean value of all theta time traces **)
meanIminmax=Mean[psdData[[All,10]]];


(*Do[
	fMax=1.0; (** found by eye where spectrum flattens out **)
	psdi=psdData\[LeftDoubleBracket]i,3,All,2\[RightDoubleBracket];
	{sFit,\[Theta]Fit,iFit,\[Theta]FitTheory,iFitTheory,sWidth,\[Theta]Width,iWidth,meanThetaFit}=calcFits[fMax,freqs,psdData\[LeftDoubleBracket]i,1,All,2\[RightDoubleBracket],psdData\[LeftDoubleBracket]i,2,All,2\[RightDoubleBracket],psdi];
	plot=Show[
		ListLinePlot[{freqs,psdi}\[Transpose],PlotRange->{{-0.01,2},All},Frame->True,FrameLabel->{"frequency f (Hz)","spectrum"},AspectRatio->1],
		Plot[iFitTheory[x],{x,0,2},PlotStyle\[Rule]Directive[Green,Opacity[0.6]],PlotRange\[Rule]All],
		PlotLabel->"mean \[Theta] fit = "<>ToString[DecimalForm[meanThetaFit,{3,2}]]<>" vs. true \[Theta] "<>ToString[DecimalForm[\[Pi]-psdData\[LeftDoubleBracket]i,9\[RightDoubleBracket],{3,2}]]
	];
	Export[FileNameJoin[{pthout,"individual_psdi_"<>ToString[Round[10^6u0Dim]]<>"_"<>IntegerString[i,10,4]<>".png"}],plot];
,{i,Round@Subdivide[1,nParallelSystems,10-1]}];*)

(*Do[
	xyz=trajectories\[LeftDoubleBracket]i,tTransient;;,{1,2,3}\[RightDoubleBracket];
	speed=Norm/@(Differences[xyz]/dt)rDropletDim 10^6;
	plot=ListLinePlot[{Range[Length[speed]]dt,speed}\[Transpose],PlotRange->All,Frame->True,FrameLabel->{"time (s)","speed (\[Mu]m/s)"},AspectRatio->1];
	Export[FileNameJoin[{pthout,"individual_speeds_"<>ToString[Round[10^6u0Dim]]<>"_"<>IntegerString[i,10,4]<>".png"}],plot];
,{i,Round@Subdivide[1,nParallelSystems,10]}];*)

(** output **)
{freqs,mPSDs,mPSD\[Theta],mPSDi,stdPSD\[Theta],stdPSDi,sData,thetaData,phiData,iData,phasesS,phases\[Theta],phasesI,
meanThetas,meanMeanTheta,thetaNoise,mThetaNoiseMean,mThetaNoiseC,meanIminmax,histoTheta}
]


(* ::Subsection::Closed:: *)
(*calcPhysicalParameters*)


(* ::Text:: *)
(*TODO : merge with calcForceFactors[]*)


(* ::Input:: *)
(*calcPhysicalParameters[u0Dim,rDropletDim];*)


(* ::Subsubsection::Closed:: *)
(*test*)


(* ::Input:: *)
(*calcPhysicalParameters[10^-6,3.5*^-6];*)
(*m*)


(* ::Input:: *)
(*mDim gDim rcomDim/(rDropletDim+rEcoliShortDim)^2*)


(* ::Input:: *)
(*rcm/(rDroplet+rEcoliShort)*)


(* ::Input:: *)
(*dragTransDim*)
(*dragRotDim*)


(* ::Input:: *)
(*8\[Pi] \[Eta]  rDropletDim*)
(*8\[Pi] \[Eta]  rDropletDim^3/(rDropletDim+rEcoliShortDim)^2*)


(* ::Subsubsection:: *)
(*code*)


calcPhysicalParameters[u0Dim_:20.0*^-6,rDropletDim_:3.5*^-6]:=Module[{},
(*u0Dim= 20.0*^-6;                            (** constant speed: 20 micron/s **)
rDropletDim=3.5*^-6;
*)

calcGeometricalParameters[rDropletDim];

\[Eta]=8.90*^-4; (** kg/s **)


dragTransDim=6\[Pi] \[Eta]  rEcoliShortDim (*(rDropletDim+rEcoliShortDim)^2*);    (** in kg m / s **)
dragRotDim=8\[Pi] \[Eta]  rDropletDim;    (** in kg m / s **)
dragDim=dragTransDim+dragRotDim;  (** total drag factor **)
swimFactor=swimMultiplier;  (** 1.0 **)
dragSwimDim=swimFactor 6\[Pi] \[Eta] rEcoliShortDim;   (** active friction **)

{\[Sigma]1,\[Sigma]2}=0.1{1,1};                     (** translational noise intensities **)
{gDim,diffR}={9.81,3.5};              (** gravitational constant in m/s^2; rotational diffusion constant in rad^2/s **)

dragFactor=dragTransDim;

(** spatially rescale units **)
{u0,g}={u0Dim,gDim}/rDropletDim;

prefactor=m g rcm/dragFactor;
(*prefactor=0;*)
]


(* ::Subsection::Closed:: *)
(*calcStateStatistics*)


calcStateStatistics[trajectories_,nSteps_,dt_]:=Module[{rtStates,rtRatios,rTimes1,rTimes2,tRates1,tRates2,residenceData,data0,data1,meanResidenceTimes,meanTransitionRates},
	rtStates=Map[If[#<0.5,0,1]&,trajectories[[All,All,7]],{2}];   (** binarized run-tumble state variable **)
	rtRatios=N[Count[#,1]/nSteps]&/@rtStates;
	{rTimes1,rTimes2,tRates1,tRates2}=Table[
		residenceData=Sort[{First[#],Length[#]}&/@Split[rtState]];
		{data0,data1}=Flatten/@(Reap[Do[Sow[d[[2]],If[d[[1]]==0,list0,list1]],{d,residenceData}],{list0,list1}][[2]]);
		meanResidenceTimes=dt If[#=={},0.0,N[Mean[#]]]&/@{data0,data1};
		meanTransitionRates=If[#==0.0,0,1/#]&/@meanResidenceTimes;
		Flatten[{meanResidenceTimes,meanTransitionRates}]
	,{rtState,rtStates}]\[Transpose];
	meanResidenceTimes=Mean/@{rTimes1,rTimes2};
	meanTransitionRates=Mean/@{tRates1,tRates2};
	(*Print[meanStd[rtRatios]];
	Print["residence times (s): ",meanResidenceTimes];
	Print["transition rates (s^-1): ",meanTransitionRates];*)
	{rtRatios,meanResidenceTimes,meanTransitionRates}
]


(* ::Subsection::Closed:: *)
(*calcTangentVelocitiesForces*)


(* ::Subsubsection::Closed:: *)
(*slow Mathematica code*)


(* ::Code:: *)
(*(** slow manual functions **)*)
(*Table[*)
(*	xyz=trajectory[[All,{1,2,3}]];    (** position **)*)
(*	pxyz=trajectory[[All,{4,5,6}]];   (** polarization **)*)
(*	s=trajectory[[All,7]];            (** RT state variable **)*)
(*	projectionMatrices=sphereSurfaceProjector[#,{0,0,0}]&/@xyz[[;;-2]];*)
(*	velocity=Differences[xyz]/dt;(*MapThread[#1.#2&,{projectionMatrices,Differences[xyz]/dt}];*)*)
(*	swimForce=swimFactor MapThread[(1.0-#1) #2&,{s[[;;-2]],pxyz[[;;-2]]}];  (** polarization vector is already tangent to sphere **)*)
(*	gravForce=(# . {0,0,1} gravFactor)&/@projectionMatrices;*)
(*	normVelocity=Norm/@velocity;*)
(*	zVelocity=velocity[[All,3]];*)
(*	zSwimForce=swimForce[[All,3]];*)
(*	zGravForce=gravForce[[All,3]];*)
(*	zSwimGravForce=zSwimForce+zGravForce;*)
(*	Join[{*)
(*		{Mean@#,Min@#,Max@#}&@normVelocity,*)
(*		{Mean[Norm/@swimForce]},{Mean[Norm/@gravForce]},{Mean[Norm/@(swimForce+gravForce)]},*)
(*		{Mean@#,Min@#,Max@#}&@zVelocity,{Mean@#,Min@#,Max@#}&@zSwimForce,{Mean@#,Min@#,Max@#}&@zGravForce,{Mean@#,Min@#,Max@#}&@zSwimGravForce*)
(*	}]*)
(*,{trajectory,trajectories}]\[Transpose];*)


(* ::Subsubsection::Closed:: *)
(*compiled code*)


(* ::Input:: *)
(*Needs["CompiledFunctionTools`"]*)
(*CompilePrint[calcTangentVelocitiesForcesC]*)


calcTangentVelocitiesForcesC=Compile[{{trajectory,_Real,2},{dt,_Real},{swimFactor,_Real},{gravFactor,_Real}},
Module[{xyz,pxyz,s,normals,projectionMatrices,velocity,swimForce,gravForce,normVelocity,zVelocity,zSwimForce,zGravForce,zSwimGravForce,output=Table[0.0,{18}]},
xyz=trajectory[[All,{1,2,3}]];    (** position **)
pxyz=trajectory[[All,{4,5,6}]];   (** polarization **)
s=trajectory[[All,7]];             (** RT state variable **)

normals=#/Norm[#]&/@xyz[[;;-2]];
projectionMatrices=({{1,0,0},{0,1,0},{0,0,1}}-{{#[[1]]#[[1]],#[[1]] #[[2]],#[[1]]#[[3]]},{#[[2]]#[[1]],#[[2]]#[[2]],#[[2]]#[[3]]},{#[[3]]#[[1]],#[[3]]#[[2]],#[[3]]#[[3]]}})&/@xyz[[;;-2]];
velocity=Differences[xyz]/dt;(*MapThread[#1.#2&,{projectionMatrices,Differences[xyz]/dt}];*)
swimForce=swimFactor MapThread[(1.0-#1) #2&,{s[[;;-2]],pxyz[[;;-2]]}];  (** polarization vector is already tangent to sphere **)
gravForce=(# . {0,0,1} gravFactor)&/@projectionMatrices;
normVelocity=Norm/@velocity;
zVelocity=velocity[[All,3]];
zSwimForce=swimForce[[All,3]];
zGravForce=gravForce[[All,3]];
zSwimGravForce=zSwimForce+zGravForce;
output[[1;;3]]={Mean@#,Min@#,Max@#}&@normVelocity;
output[[4;;6]]={Mean[Norm/@swimForce],Mean[Norm/@gravForce],Mean[Norm/@(swimForce+gravForce)]};
output[[7;;9]]={Mean@#,Min@#,Max@#}&@zVelocity;
output[[10;;12]]={Mean@#,Min@#,Max@#}&@zSwimForce;
output[[13;;15]]={Mean@#,Min@#,Max@#}&@zGravForce;
output[[16;;18]]={Mean@#,Min@#,Max@#}&@zSwimGravForce;
output
],CompilationTarget->"C",RuntimeOptions->"Speed",RuntimeAttributes->{Listable},Parallelization->True];


(* ::Subsubsection::Closed:: *)
(*code*)


calcTangentVelocitiesForces[trajectories_,dt_,swimFactor_,gravFactor_,rDropletDim_]:=Module[{velocityForceData,xyz,pxyz,s,projectionMatrices,velocity,swimForce,gravForce,normVelocity,zVelocity,zSwimForce,zGravForce,zSwimGravForce},

velocityForceData=calcTangentVelocitiesForcesC[trajectories,dt,swimFactor,gravFactor];

(** multiply velocity (m/s) and force (kg m/s^2) by space rescaling factor to go from numerical to physical values **)
velocityForceData*=rDropletDim;

Flatten[meanStd/@velocityForceData]
]


(* ::Subsection::Closed:: *)
(*createMovieFramesOnlyBacterium*)


(* ::Text:: *)
(*requires droplet_optics code to be run in advance to compute liquid crystal alignments*)


(* ::Subsubsection::Closed:: *)
(*debug*)


(* ::Input:: *)
(*origin={0.8 1.0,0.8*-1.0,-1};*)
(*plotAxesCross=Graphics3D[{Black,Arrowheads[0.03],Arrow[Tube[{origin,origin+0.35#}]]&/@{{-1,0,0},{0,1,0},{0,0,0.9 1}}}];*)
(*Show[*)
(*Graphics3D[{Opacity[0.1],Sphere[]}],*)
(*plotAxesCross*)
(*,Axes->False*)
(*,Boxed->False*)
(*,ImageSize->500*)
(*,PlotRange->1.4{{-1,1},{-1,1},{-1,1}}*)
(*,ViewVertical->{0.376,0.248,0.893}*)
(*,ViewPoint->{2.813,1.728,0.742}*)
(*]*)


(* ::Subsubsection:: *)
(*code*)


createMovieFramesOnlyBacterium[trajectory_,i_,dt_]:=Module[{color,counter,outputDirectory,rtData,position,theta,phi,
orientation,front,end,runColor,tumbleColor,plotBacteriumSphere,bacteriumAge,timeLabel,plot,liquidCrystals,droplet,plotDropletLC},

color={colorYoung,colorOld}[[i]];
{runColor,tumbleColor}={color,Darker[color,0.5]};

counter=0;
outputDirectory=FileNameJoin[{NotebookDirectory[],"mmamovie_"<>ToString[i]}];
If[!DirectoryQ[#],CreateDirectory[#]]&@outputDirectory;


positions=trajectory[[All,1;;3]];
thetas=ArcCos/@positions[[All,3]];   (** 0 - top pole (z=+1), \[Pi] - bottom pole (z=-1) **)
phis=phaseUnwrap[ArcTan[#[[1]],#[[2]]]&/@positions[[All,{1,2}]]];
orientations=trajectory[[All,4;;6]];
rtDatas=trajectory[[All,7]];


Do[
	position=positions[[t]];
	orientation=orientations[[t]];
	theta=thetas[[t]];
	phi=phis[[t]];
	rtData=If[#<0.5,0,1]&/@rtDatas[[;;t]];
	{front,end}=Map[(1+rEcoliShort)position+# ((*rEcoliLong*)2rEcoliShort-rEcoliShort) orientation&,{1,-1}];
	
	(** output plot **)
	liquidCrystals=renderLCs[solution,{theta,phi}];

	droplet=ParametricPlot3D[{
			Evaluate[RotationMatrix[{{0,0,1},position}] . er[\[Theta],\[Phi]]]
			,Evaluate[RotationMatrix[{{0,0,1},position}] . er[\[Pi]-\[Theta],\[Phi]]]
		},{\[Theta],0,0.5\[Pi]},{\[Phi],0,2\[Pi]}
		,Mesh->None
		,Axes->False
		,Boxed->False
		,PlotPoints->50
		,PlotStyle->Evaluate[Directive[Opacity[0.2],Specularity[White,65],#]&/@{Lighter[Red],Lighter[Gray,0.7]}]
	];

	plotDropletLC=Show[droplet,liquidCrystals,PlotRange->All,Lighting->"Neutral"];
	(*Graphics3D[{Opacity[0.1],Sphere[]}]*)
	
	(** at origin **)
	(*plotAxesCross=Graphics3D[{Black,Arrow[Tube[{{0,0,0},1.3#}]]&/@{{1,0,0},{0,1,0},{0,0,1}}}];*)
	
	(** bottom left **)
	origin={0.8 1.0,0.8*-1.0,-1};
	plotAxesCross=Graphics3D[{Black,Arrowheads[0.03],Arrow[Tube[{origin,origin+0.35#}]]&/@{{-1,0,0},{0,1,0},{0,0,0.9 1}}}];

	plotBacteriumSphere=Show[
		plotDropletLC
		,plotAxesCross
		,Graphics3D[{ColorData[97,1],AbsoluteThickness[5],Line[#&/@trajectory[[;;t,1;;3]],VertexColors->(If[#==0,runColor,tumbleColor]&/@rtData)]}]	
		,Graphics3D[{If[rtData[[t]]==0,runColor,tumbleColor],CapsuleShape[{front,end},rEcoliShort]}]
		,Axes->False
		,Boxed->False
		,ImageSize->500
		,PlotRange->1.4{{-1,1},{-1,1},{-1,1}}
		,ViewVertical->{0.376,0.248,0.893}
		,ViewPoint->{2.813,1.728,0.742}
	];

	bacteriumAge=If[i==1,"Fresh E. coli","Starved E. coli"];
	timeLabel=Style[bacteriumAge<>"\ntime \!\(\*
	StyleBox[\"t\",\nFontSlant->\"Italic\"]\) = "<>ToString@DecimalForm[t dt,{5,Round[-Log10[dt]]}]<>"\[ThinSpace]s",FontFamily->"Helvetica",26];
	
	(** output plot **)
	plot=Grid[{
		{timeLabel},
		{plotBacteriumSphere}
	}];
	
	Export[FileNameJoin[{outputDirectory,"p_"<>IntegerString[counter++,10,IntegerLength[Length[trajectory]]+1]<>".png"}],plot];
,{t,1,Length[trajectory],1}]
]


(* ::Subsection::Closed:: *)
(*createMovieFramesPositionOrientation*)


(* ::Text:: *)
(*position + orientation spheres*)


createMovieFramesPositionOrientation[trajectory_,i_,dt_]:=Module[{color,counter,outputDirectory,rtData,position,orientation,front,end,runColor,tumbleColor,bacteriumAge,timeLabel,plotBacteriumSphere,plot},

color=ColorData[97,i];
counter=0;
outputDirectory=FileNameJoin[{NotebookDirectory[],"mmamovie_"<>ToString[i]}];
If[!DirectoryQ@#,CreateDirectory@#]&@outputDirectory;
Do[

rtData=If[#<0.5,0,1]&/@trajectory[[;;t,7]];
position=trajectory[[t,1;;3]];
orientation=trajectory[[t,4;;6]];
{front,end}=Map[(1+rEcoliShort)position+# (rEcoliLong-rEcoliShort) orientation&,{1,-1}];
{runColor,tumbleColor}={color,Darker[color,0.5]};

plotBacteriumSphere=Show[
Graphics3D[{Opacity[0.1],Sphere[]}],
Graphics3D[{ColorData[97,1],AbsoluteThickness[3],Line[#&/@trajectory[[;;t,1;;3]],VertexColors->(If[#==0,runColor,tumbleColor]&/@rtData)]}],	
Graphics3D[{If[rtData[[t]]==0,runColor,tumbleColor],CapsuleShape[{front,end},rEcoliShort]}],
Graphics3D[{Black,Arrow[Tube[1.3{{0,0,0},#}]]&/@{{1,0,0},{0,1,0},{0,0,1}}}]
,Boxed->False,Axes->False,ImageSize->500,PlotRange->1.4{{-1,1},{-1,1},{-1,1}}
,ViewVertical->{0.376,0.248,0.893},ViewPoint->{2.813,1.728,0.742}];

bacteriumAge=If[i==1,"Young E. coli","Aged E. coli"];
timeLabel=Style[bacteriumAge<>"\ntime t = "<>ToString@DecimalForm[t dt,{5,Round[-Log10[dt]]}],FontFamily->"Helvetica",26];

(** output plot **)
plot=Grid[{
{timeLabel},
{plotBacteriumSphere}
}];


Export[FileNameJoin[{outputDirectory,"p_"<>IntegerString[counter++,10,IntegerLength[Length[trajectory]]+1]<>".png"}],plot]
,{t,1,Length[trajectory],1}]
]


(* ::Subsection::Closed:: *)
(*doNoiseAnalysisPlots*)


doNoiseAnalysisPlots[plotNoiseQ_,thetaNoise_,mThetaNoiseC_,thetaData_]:=Module[{plotOps,plotNoise,plotNoiseCorrelation,plotThetaNoise},
If[plotNoiseQ,
plotOps={AspectRatio->1,Frame->True,PlotRange->All,ImageSize->{Automatic,400}};
plotNoise=ListPlot[thetaNoise,Evaluate[plotOps],FrameLabel->{"time t","\[Theta] noise"},Joined->True];
plotNoiseCorrelation=ListPlot[mThetaNoiseC,Evaluate[plotOps],FrameLabel->{"time t","\[Theta] noise autocorrelation"},Joined->True];
plotThetaNoise=ListPlot[{thetaData[[;;-2]],thetaNoise}\[Transpose],Evaluate[plotOps],FrameLabel->{"\[Theta]","\[Theta] noise"}];
Print[Grid[{{plotNoise,plotNoiseCorrelation,plotThetaNoise}}]];
]
]


(* ::Subsection::Closed:: *)
(*doOverviewPlotsSaveData*)


(* ::Input:: *)
(*doOverviewPlotsSaveData[rDropletDim,simData,pthout,u0Dims,spectraTheta,spectraI,freqs,fMax,u0Max];*)


(* ::Subsubsection:: *)
(*debug*)


(* ::Subsubsection:: *)
(*code*)


doOverviewPlotsSaveData[parameterString_,simData_,pthout_,u0Dims_,spectraTheta_,spectraI_,freqs_,fMax_,u0Max_]:=Module[
{plotOps,plotWidth,plotResidenceTimes,plotMeanTheta,allSpectraMax,plotSpectrumTheta,
plotSpectrumI,plotSpectrumThetaColorbar,plotSpectrumIColorbar,fMaxIndex,fTicks,plotOpsSpectra,u0DimData,iWidthData,
meanRunTimeData,meanTumbleTimeData,meanThetaData,meanThetaFitData,fMatrix,u0Matrix,data,u0},

(** overview plots **)
plotOps={AspectRatio->1,Frame->True,FrameStyle->Directive[Black,28,AbsoluteThickness[2]],ImageSize->500,
PlotRangePadding->{Scaled[0.02],{0,Scaled@0.05}},ColorFunction->"BlueGreenYellow",Mesh->All,ScalingFunctions->{"Reverse",Identity}};

u0DimData=10^6 GeneralUtilities`AssociationTranspose[simData]["u0Dim"];
iWidthData=GeneralUtilities`AssociationTranspose[simData]["iWidth"];
plotWidth=ListLinePlot[{u0DimData,iWidthData}\[Transpose],Evaluate[plotOps],PlotRange->{10^6 MinMax[u0Dims],{0,All}},FrameLabel->{"speed \!\(\*
StyleBox[SubscriptBox[\"u\", \"0\"],\nFontSlant->\"Italic\"]\) (\!\(\*SuperscriptBox[\(10\), \(-6\)]\) m/s)","width \!\(\*
StyleBox[\"w\",\nFontSlant->\"Italic\"]\) (Hz)"}];
Export[FileNameJoin[{pthout,"overview_width"<>parameterString<>"."<>#}],plotWidth]&/@{"png"};
Print[plotWidth];

meanRunTimeData=GeneralUtilities`AssociationTranspose[simData]["meanRunTime"];
meanTumbleTimeData=GeneralUtilities`AssociationTranspose[simData]["meanTumbleTime"];
plotResidenceTimes=ListPlot[{{u0DimData,meanRunTimeData}\[Transpose],{u0DimData,meanTumbleTimeData}\[Transpose]},
Evaluate[plotOps],PlotRange->{{1,u0Max},{0,1.0}},FrameLabel->{"speed \!\(\*
StyleBox[SubscriptBox[\"u\", \"0\"],\nFontSlant->\"Italic\"]\) (\!\(\*SuperscriptBox[\(10\), \(-6\)]\) m/s)","residence times (s)"}];
Export[FileNameJoin[{pthout,"overview_ratios"<>parameterString<>"."<>#}],plotResidenceTimes]&/@{"png"};
(*Print[Grid[{{plotWidth,plotResidenceTimes}},Spacings->{1, 1}]];*)

(** polar angle plots **)
fitMeanThetaQ=False;
meanThetaData=GeneralUtilities`AssociationTranspose[simData]["meanTheta"];
If[fitMeanThetaQ,
	meanThetaFitData=\[Pi]-GeneralUtilities`AssociationTranspose[simData]["meanThetaFit"];
	plotMeanTheta=ListLinePlot[{{u0DimData,meanThetaData}\[Transpose],{u0DimData,meanThetaFitData}\[Transpose]},ColorFunction->Automatic,ScalingFunctions->{"Reverse",Identity},Evaluate[plotOps],
		PlotRange->{{1,u0Max},{0,\[Pi]}},FrameLabel->{"speed \!\(\*
StyleBox[SubscriptBox[\"u\", \"0\"],\nFontSlant->\"Italic\"]\) (\!\(\*SuperscriptBox[\(10\), \(-6\)]\) m/s)","mean angle \!\(\*
StyleBox[\"\[LeftAngleBracket]\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"\[Theta]\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"\[RightAngleBracket]\",\nFontSlant->\"Italic\"]\) (rad)"},PlotLegends->{"true \[LeftAngleBracket]\[Theta]\[RightAngleBracket]","fitted \[LeftAngleBracket]\[Theta]\[RightAngleBracket]"}
	];
,
	plotMeanTheta=Show[
		ListLinePlot[{u0DimData,180.0/\[Pi] meanThetaData}\[Transpose],ColorFunction->Automatic,ScalingFunctions->{"Reverse",Identity},Evaluate[plotOps],
			PlotRange->{10^6 MinMax[u0Dims],{0,All}},FrameLabel->{"speed \!\(\*
StyleBox[SubscriptBox[\"u\", \"0\"],\nFontSlant->\"Italic\"]\) (\!\(\*SuperscriptBox[\(10\), \(-6\)]\) m/s)","mean angle \!\(\*
StyleBox[\"\[LeftAngleBracket]\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"\[Theta]\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"\[RightAngleBracket]\",\nFontSlant->\"Italic\"]\) (\[Degree])"},PlotLegends->{"sim. \[LeftAngleBracket]\[Theta]\[RightAngleBracket]"}
		],
		Plot[180.0/\[Pi] maxTheta[10^-6 u0],{u0,0,10^6 Max[u0Dims]},ScalingFunctions->{"Reverse",Identity},PlotStyle->Directive[Dashed,Black],PlotLegends->{"theo. \!\(\*SubscriptBox[\(\[Theta]\), \(max\)]\)"}]
	,PlotRange->All];
];
Export[FileNameJoin[{pthout,"overview_meanTheta_cVals"<>parameterString<>"."<>#}],plotMeanTheta]&/@{"png"};
Print[plotMeanTheta];


(** spectra plots **)
fMaxIndex=FirstPosition[freqs,_?(#>=fMax&)][[1]];
fTicks=Table[i{20,0.1},{i,1,5}];
(*allSpectraMax=Max[spectraTheta\[LeftDoubleBracket]All,;;fMaxIndex\[RightDoubleBracket]\[Transpose]];*)
plotOpsSpectra={FrameStyle->Directive[Black,20],InterpolationOrder->0,ColorFunction->"BlueGreenYellow",
FrameLabel->{"speed \!\(\*
StyleBox[SubscriptBox[\"u\", \"0\"],\nFontSlant->\"Italic\"]\) (\!\(\*SuperscriptBox[\(10\), \(-6\)]\) m/s)","frequency \!\(\*
StyleBox[\"f\",\nFontSlant->\"Italic\"]\) (Hz)"},PlotRange->All,PlotRangePadding->None,ScalingFunctions->{"Reverse",Identity},ImageSize->{Automatic,500},
ImagePadding->{{70,5},{60,10}}};

fMatrix=ConstantArray[freqs[[;;fMaxIndex]],Length[u0Dims]]\[Transpose];
u0Matrix=ConstantArray[10^6 u0Dims,fMaxIndex];
data=spectraTheta[[All,;;fMaxIndex]]\[Transpose];
plotSpectrumTheta=Show[
	ListDensityPlot[Flatten[Transpose[{u0Matrix,fMatrix,data},{3,2,1}],1],Evaluate[plotOpsSpectra]],
	Plot[avgFreqThetaNeuron[u0 10^-6],{u0,10^6 Min[u0Dims],10^6 Max[u0Dims]},PlotStyle->Red,ScalingFunctions->{"Reverse",Identity}]
];
plotSpectrumThetaColorbar=Row[{
	plotSpectrumTheta,
	Show[colorbar[{0,1},{{0,0},{0.5,""},{1,1}},False,20,500,"PSD",2,"BlueGreenYellow"],ImagePadding->{{5,50},{60,10}}]
}];


data=spectraI[[All,;;fMaxIndex]]\[Transpose];
plotSpectrumI=Show[
	ListDensityPlot[Flatten[Transpose[{u0Matrix,fMatrix,data},{3,2,1}],1],Evaluate[plotOpsSpectra]],
	Plot[{1,2}avgFreqThetaNeuron[u0 10^-6],{u0,10^6 Min[u0Dims],10^6 Max[u0Dims]},PlotStyle->{Directive[Red,AbsoluteThickness[2]],Directive[Red,Dashed,AbsoluteThickness[2]]},ScalingFunctions->{"Reverse",Identity}]
];
plotSpectrumIColorbar=Row[{
	plotSpectrumI,
	Show[colorbar[{0,1},{{0,0},{0.5,""},{1,1}},False,20,500,"PSD",2,"BlueGreenYellow"],ImagePadding->{{5,50},{60,10}}]
}];

(** plots **)
Print[Grid[{{plotSpectrumThetaColorbar,plotSpectrumIColorbar}}]];
Export[FileNameJoin[{pthout,"spectrumTheta"<>parameterString<>"."<>#}],plotSpectrumThetaColorbar]&/@{"png"};
Export[FileNameJoin[{pthout,"spectrumIntensity"<>parameterString<>"."<>#}],plotSpectrumIColorbar]&/@{"png"};

(** save data **)
Export[FileNameJoin[{pthout,"spectraIntensity"<>parameterString<>".h5"}],spectraI];
Export[FileNameJoin[{pthout,"spectraTheta"<>parameterString<>".h5"}],spectraTheta];
Export[FileNameJoin[{pthout,"psdFreqs"<>parameterString<>".h5"}],freqs];

Export[FileNameJoin[{pthout,"simAnalysisData"<>parameterString<>".h5"}],simData];
]


(* ::Subsection::Closed:: *)
(*doOverviewPlotsSaveDataRT*)


doOverviewPlotsSaveDataRT[parameterString_,simData_,pthout_,u0Dims_,spectraTheta_,spectraI_,freqs_,fMax_,u0Max_]:=Module[
{plotOps,plot\[Theta]Width,plotiWidth,plotResidenceTimes,plotMeanTheta,allSpectraMax,plotSpectrumTheta,
plotSpectrumI,plotSpectrumThetaColorbar,plotSpectrumIColorbar,fMaxIndex,fTicks,plotOpsSpectra,u0DimData,\[Theta]WidthData,iWidthData,
runtimeData,tumbletimeData,meanRunTimeData,meanTumbleTimeData,meanThetaData,meanThetaFitData,fMatrix,u0Matrix,data},

{widthMin,widthMax}={0.0,0.35};
{thetaMin,thetaMax}={0.0,45};

(** overview plots **)
plotOps={PlotRange->All,InterpolationOrder->0,FrameLabel->{"runtime \!\(\*
StyleBox[SubscriptBox[\"t\", \"R\"],\nFontSlant->\"Italic\"]\) (s)","tumbletime \!\(\*
StyleBox[SubscriptBox[\"t\", \"T\"],\nFontSlant->\"Italic\"]\) (s)"},PlotRangePadding->None,FrameStyle->Directive[Black,20,AbsoluteThickness[2]],ScalingFunctions->{"Log10","Log10"}};

runtimeData=GeneralUtilities`AssociationTranspose[simData]["tRun"];
tumbletimeData=GeneralUtilities`AssociationTranspose[simData]["tTumble"];

If[Length[runtimeData]>1\[And]Length[tumbletimeData]>1,
	
	\[Theta]WidthData=GeneralUtilities`AssociationTranspose[simData]["thetaWidth"];
	plot\[Theta]Width=ListDensityPlot[{runtimeData,tumbletimeData,\[Theta]WidthData}\[Transpose],Evaluate[plotOps],
	ColorFunctionScaling->False,ColorFunction->(GrayLevel[Rescale[#,{widthMin,widthMax}]]&),
	PlotLegends->BarLegend[Automatic,LegendLabel->"width \!\(\*SubscriptBox[\(PSD\), \(\[Theta]\)]\) (Hz)"]
	];
	Export[FileNameJoin[{pthout,"overview_thetawidth"<>parameterString<>"."<>#}],plot\[Theta]Width]&/@{"png"};
	Print[plot\[Theta]Width];

	iWidthData=GeneralUtilities`AssociationTranspose[simData]["iWidth"];
	plotiWidth=ListDensityPlot[{runtimeData,tumbletimeData,iWidthData}\[Transpose],Evaluate[plotOps],
	ColorFunctionScaling->False,ColorFunction->(GrayLevel[Rescale[#,{widthMin,widthMax}]]&),
	PlotLegends->BarLegend[Automatic,LegendLabel->"width \!\(\*SubscriptBox[\(PSD\), \(I\)]\) (Hz)"]
	];
	Export[FileNameJoin[{pthout,"overview_iwidth"<>parameterString<>"."<>#}],plotiWidth]&/@{"png"};
	Print[plotiWidth];
	
	iminmaxData=GeneralUtilities`AssociationTranspose[simData]["iminmax"];
	plotMinMax=ListDensityPlot[{runtimeData,tumbletimeData,iminmaxData}\[Transpose],Evaluate[plotOps],
	(*ColorFunctionScaling->False,ColorFunction->(GrayLevel[Rescale[#,{widthMin,widthMax}]]&),*)
	PlotLegends->BarLegend[Automatic,LegendLabel->"variation \[CapitalDelta]I"]
	];
	Export[FileNameJoin[{pthout,"overview_iminmax"<>parameterString<>"."<>#}],plotMinMax]&/@{"png"};
	
	meanRunTimeData=GeneralUtilities`AssociationTranspose[simData]["meanRunTime"];
	meanTumbleTimeData=GeneralUtilities`AssociationTranspose[simData]["meanTumbleTime"];
	plotTestRuntime=ListDensityPlot[{runtimeData,tumbletimeData,meanRunTimeData}\[Transpose],PlotRange->All,InterpolationOrder->0,ColorFunction->GrayLevel,FrameLabel->{"runtime \!\(\*
StyleBox[SubscriptBox[\"t\", \"R\"],\nFontSlant->\"Italic\"]\) (s)","tumbletime \!\(\*
StyleBox[SubscriptBox[\"t\", \"T\"],\nFontSlant->\"Italic\"]\) (s)"},PlotLegends->BarLegend[Automatic,LegendLabel->"run time \[LeftAngleBracket]\!\(\*SubscriptBox[\(t\), \(r\)]\)\[RightAngleBracket] (s)"],PlotRangePadding->None,FrameStyle->Directive[Black,20,AbsoluteThickness[2]],ScalingFunctions->{"Log10","Log10"}];
	plotTestTumbletime=ListDensityPlot[{runtimeData,tumbletimeData,meanTumbleTimeData}\[Transpose],PlotRange->All,InterpolationOrder->0,ColorFunction->GrayLevel,FrameLabel->{"runtime \!\(\*
StyleBox[SubscriptBox[\"t\", \"R\"],\nFontSlant->\"Italic\"]\) (s)","tumbletime \!\(\*
StyleBox[SubscriptBox[\"t\", \"T\"],\nFontSlant->\"Italic\"]\) (s)"},PlotLegends->BarLegend[Automatic,LegendLabel->"tumble time \[LeftAngleBracket]\!\(\*SubscriptBox[\(t\), \(t\)]\)\[RightAngleBracket] (s)"],PlotRangePadding->None,FrameStyle->Directive[Black,20,AbsoluteThickness[2]],ScalingFunctions->{"Log10","Log10"}];
	Export[FileNameJoin[{pthout,"overview_test_runtime"<>parameterString<>"."<>#}],plotTestRuntime]&/@{"png"};
	Export[FileNameJoin[{pthout,"overview_test_tumbletime"<>parameterString<>"."<>#}],plotTestTumbletime]&/@{"png"};

	meanThetaData=GeneralUtilities`AssociationTranspose[simData]["meanTheta"];
	meanThetaFitData=\[Pi]-GeneralUtilities`AssociationTranspose[simData]["meanThetaFit"];
	plotMeanTheta=ListDensityPlot[{runtimeData,tumbletimeData,180/\[Pi] meanThetaData}\[Transpose],PlotRange->All,InterpolationOrder->0,ColorFunctionScaling->False,ColorFunction->(GrayLevel[Rescale[#,{thetaMin,thetaMax}]]&),FrameLabel->{"runtime \!\(\*
StyleBox[SubscriptBox[\"t\", \"R\"],\nFontSlant->\"Italic\"]\) (s)","tumbletime \!\(\*
StyleBox[SubscriptBox[\"t\", \"T\"],\nFontSlant->\"Italic\"]\) (s)"},PlotLegends->BarLegend[Automatic,LegendLabel->"mean angle \!\(\*
StyleBox[\"\[LeftAngleBracket]\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"\[Theta]\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"\[RightAngleBracket]\",\nFontSlant->\"Italic\"]\) (\[Degree])"],PlotRangePadding->None,FrameStyle->Directive[Black,20,AbsoluteThickness[2]],ScalingFunctions->{"Log10","Log10"}];
	plotMeanThetaFit=ListDensityPlot[{runtimeData,tumbletimeData,180/\[Pi] meanThetaFitData}\[Transpose],PlotRange->All,InterpolationOrder->0,ColorFunction->GrayLevel,FrameLabel->{"runtime \!\(\*
StyleBox[SubscriptBox[\"t\", \"R\"],\nFontSlant->\"Italic\"]\) (s)","tumbletime \!\(\*
StyleBox[SubscriptBox[\"t\", \"T\"],\nFontSlant->\"Italic\"]\) (s)"},PlotLegends->BarLegend[Automatic,LegendLabel->"mean angle fit \!\(\*
StyleBox[\"\[LeftAngleBracket]\",\nFontSlant->\"Italic\"]\)\!\(\*
StyleBox[\"\[Theta]\",\nFontSlant->\"Italic\"]\)\!\(\*SubscriptBox[
StyleBox[\"\[RightAngleBracket]\",\nFontSlant->\"Italic\"], \(fit\)]\) (\[Degree])"],PlotRangePadding->None,FrameStyle->Directive[Black,20,AbsoluteThickness[2]],ScalingFunctions->{"Log10","Log10"}];
	Export[FileNameJoin[{pthout,"overview_meanTheta"<>parameterString<>"."<>#}],plotMeanTheta]&/@{"png"};
	Export[FileNameJoin[{pthout,"overview_meanThetaFit"<>parameterString<>"."<>#}],plotMeanThetaFit]&/@{"png"};
	Print[plotMeanTheta];

	(** save data **)
	Export[FileNameJoin[{pthout,"simAnalysisData"<>parameterString<>".h5"}],simData];
];

]


(* ::Subsection::Closed:: *)
(*doPlots*)


(* ::Subsubsection::Closed:: *)
(*code*)


Clear[doPlots];
doPlots[parameterString_,times_,trajectory_,tfinalSim_,freqs_,mPSDs_,sFit_,\[Theta]Fit_,iFit_,\[Theta]FitTheory_,iFitTheory_,
mPSD\[Theta]_,mPSDi_,stdPSD\[Theta]_,stdPSDi_,meanThetas_,rEcoli_,pthout_,observables_,debugPlotQ_,histoTheta_,fMax_:1.0]:=Module[
{sData,thetaData,phiData,iData,ops,pPSDs,plot,sPlot,\[Theta]Plot,\[Phi]Plot,iPlot,pPSD\[Theta],pPSDi,plotSpectra,runColor,tumbleColor,rtData,position,orientation,tPlot,front,end,
descriptions,numericalStrings0,numericalStrings,plotData,tMax,plotTrajectory,meanThetasPlot,debugQ,thetaDataWrapped,phiDataWrapped,plotPiecesQ,histoThetaIComparisonPlot},

\[Theta]MaxDeg=(180.0/\[Pi] maxTheta[u0Dim,fudgeFactor]);

(** debug **)
plotThetaIComparisonQ=False;
If[plotThetaIComparisonQ,
	histoThetaIComparisonPlot=Show[
		ListLinePlot[{Subdivide[0.0,180,Length[#]-1],#/Max[#]}\[Transpose]&@histoTheta,InterpolationOrder->0
		,AspectRatio->1,ImageSize->{Automatic,400}
		,Frame->True,FrameStyle->Directive[Black,22,AbsoluteThickness[2]],FrameLabel->{"polar angle \!\(\*
	StyleBox[\"\[Theta]\",\nFontSlant->\"Italic\"]\)\!\(\*
	StyleBox[\" \",\nFontSlant->\"Italic\"]\)(\[Degree])","occurrence"}
		,Filling->Bottom,PlotRangePadding->{{0,0},{0,Scaled[0.05]}},PlotRange->All
		,If[Im[\[Theta]MaxDeg]==0,Epilog->{Red,Dashed,InfiniteLine[{{#,0},{#,1}}]&@\[Theta]MaxDeg},Sequence@@{}]],
		Plot[thetaPhiIntensityInterpolation[-Abs[\[Pi]/180\[Theta]-\[Pi]/2.0]+\[Pi]/2.0,0],{\[Theta],0,180},PlotRange->All,PlotStyle->ColorData[97,2]]
	];
	Export[FileNameJoin[{pthout,"histoThetaIComparison"<>parameterString<>".png"}],histoThetaIComparisonPlot];
];

(** debug **)
plotMeanThetasQ=False;
If[plotMeanThetasQ,
	meanThetasPlot=Histogram[meanThetas,{0,\[Pi],0.05},PlotRange->{{0,\[Pi]},All},AspectRatio->1,ImageSize->{Automatic,400},Frame->True,FrameLabel->{"mean angle \[LeftAngleBracket]\[Theta]\[RightAngleBracket]","occurrence"}];
	Export[FileNameJoin[{pthout,"meanThetas"<>parameterString<>".png"}],meanThetasPlot];
];

sData=trajectory[[All,7]];
thetaData=ArcCos/@trajectory[[All,3]];
phiData=phaseUnwrap[ArcTan[#[[1]],#[[2]]]&/@trajectory[[All,{1,2}]]];
(*iData=Cos[#]^2&/@thetaData;*)
(** interpolation mirrored around \[Pi]/2 and \[Pi]/4 for \[Theta] and \[Phi] **)
thetaDataWrapped=-Abs[thetaData-\[Pi]/2.0]+\[Pi]/2.0;
phiDataWrapped=-Abs[Mod[phiData,\[Pi]/2.0]-\[Pi]/4.0]+\[Pi]/4.0;
iData=MapThread[thetaPhiIntensityInterpolation[#1,#2]&,{thetaDataWrapped,phiDataWrapped}];


ops={AspectRatio->1,ImageSize->400,Frame->True,FrameStyle->Directive[Black,20,AbsoluteThickness[2]],ColorFunctionScaling->False,PlotRangePadding->{0,Scaled[0.05]},ImagePadding->{{70,15},{60,10}}};
sPlot=ListLinePlot[{times,sData}\[Transpose],Evaluate[ops],FrameLabel->{"Time \!\(\*
StyleBox[\"t\",\nFontSlant->\"Italic\"]\) (s)","RT state \!\(\*
StyleBox[\"s\",\nFontSlant->\"Italic\"]\)"},PlotRange->{{0,60},{0,1}}];
\[Theta]Plot=ListLinePlot[{times,thetaData}\[Transpose],Evaluate[ops],FrameLabel->{"Time \!\(\*
StyleBox[\"t\",\nFontSlant->\"Italic\"]\) (s)","\!\(\*
StyleBox[\"\[Theta]\",\nFontSlant->\"Italic\"]\)"},Prolog->{Gray,Line[{{0,#},{tfinalSim,#}}]&@(0.5\[Pi])},PlotRange->{{0,60},(*{\[Pi],0}*){0,\[Pi]}},ScalingFunctions->{Identity,Identity(*"Reverse"*)}
,Prolog->{Gray,Dashed,Line[{{0.0,#},{tfinalSim,#}}]&@maxTheta[u0Dim]}];
\[Phi]Plot=ListLinePlot[{times,phiData}\[Transpose],Evaluate[ops],FrameLabel->{"Time \!\(\*
StyleBox[\"t\",\nFontSlant->\"Italic\"]\) (s)","\!\(\*
StyleBox[\"\[Phi]\",\nFontSlant->\"Italic\"]\)"},PlotRange->{{0,60},All}];
iPlot=ListLinePlot[{times,iData}\[Transpose],Evaluate[ops],FrameLabel->{"Time \!\(\*
StyleBox[\"t\",\nFontSlant->\"Italic\"]\) (s)","Intensity \!\(\*
StyleBox[\"I\",\nFontSlant->\"Italic\"]\)"},PlotRange->{{0,60},{0,All}(*{0,1}*)}];

pPSDs=Show[
	ListLinePlot[{freqs,mPSDs}\[Transpose],Evaluate[ops],FrameLabel->{"Frequency \!\(\*
StyleBox[\"f\",\nFontSlant->\"Italic\"]\) (Hz)","PSD(s)"},PlotRange->{{0,fMax},All}],
	Plot[sFit[x],{x,0,fMax},PlotStyle->Red,PlotRange->All]
	(*,Plot[2(1)^2/(Total[meanResidenceTimes](Total[meanTransitionRates]^2+(2\[Pi] f)^2)),{f,0,fMax},PlotStyle\[Rule]Darker@Green,PlotRange\[Rule]All]*)
	(** Eq. 1.49 in Lindner's chapter for telegraph noise ~ to double potential noise spectrum / amplitudes are wrong, since non-Markovian effects are neglected **)
];
pPSD\[Theta]=Show[
	ListLinePlot[{{freqs,mPSD\[Theta]}\[Transpose],{freqs,mPSD\[Theta]-stdPSD\[Theta]}\[Transpose],{freqs,mPSD\[Theta]+stdPSD\[Theta]}\[Transpose]},Evaluate[ops],FrameLabel->{"Frequency \!\(\*
StyleBox[\"f\",\nFontSlant->\"Italic\"]\) (Hz)","PSD(\[Theta])"},PlotRange->{{0,fMax},All},PlotStyle->(Directive[ColorData[97,1],AbsoluteThickness[#]]&/@{4,1,1}),Filling->{3->{2}}]
	(*,Plot[\[Theta]Fit[x],{x,0,fMax},PlotStyle\[Rule]Directive[Red,Opacity[0.6]],PlotRange\[Rule]All]*)
	,Plot[Max[mPSD\[Theta]]\[Theta]FitTheory[f],{f,0,fMax},PlotStyle->Directive[Green,Opacity[0.6]],PlotRange->All]
	,PlotLegends->{"data","\!\(\*SuperscriptBox[\(Lorentz\), \(2\)]\) fit","theory fit"}
];
pPSDi=Show[
	ListLinePlot[{{freqs,mPSDi}\[Transpose],{freqs,mPSDi-stdPSDi}\[Transpose],{freqs,mPSDi+stdPSDi}\[Transpose]},Evaluate[ops],FrameLabel->{"Frequency \!\(\*
StyleBox[\"f\",\nFontSlant->\"Italic\"]\) (Hz)","PSD(i)"},PlotRange->{{0,fMax},All},PlotStyle->(Directive[ColorData[97,1],AbsoluteThickness[#]]&/@{4,1,1}),Filling->{3->{2}}]
	(*,Plot[iFit[x],{x,0,fMax},PlotStyle\[Rule]Directive[Red,Opacity[0.6]],PlotRange\[Rule]All]*)
	,Plot[Max[mPSDi]iFitTheory[f],{f,0,fMax},PlotStyle->Directive[Green,Opacity[0.6]],PlotRange->All]
	,PlotLegends->{"data","\!\(\*SuperscriptBox[\(Lorentz\), \(2\)]\) fit","theory fit"}
];

plotPiecesQ=False;
plotSpectra=Grid[{(*{sPlot,pPSDs},*){\[Theta]Plot,pPSD\[Theta]},{iPlot,pPSDi}}];
If[plotPiecesQ\[And]pthout!="",Export[FileNameJoin[{pthout,"timeseries"<>parameterString<>"."<>#}],plotSpectra]&/@{"png"}];

descriptions={"width(\[Theta])","width(i)","mean \!\(\*SubscriptBox[\(T\), \(t\)]\)/\!\(\*SubscriptBox[\(T\), \(all\)]\) ratio","mean res time R","mean res time T","mean res rate R","mean res rate T"};
numericalStrings0=ToString[DecimalForm[#,{4,3}]]&/@observables;
numericalStrings=Flatten@Join[#[[1;;2]],{#[[3]]<>"\[PlusMinus]"<>#[[4]]},#[[5;;]]]&@numericalStrings0;
plotData=Grid[{descriptions,numericalStrings}\[Transpose],Frame->All];

(** trajectory on sphere **)
sData=trajectory[[All,7]];
rtData=If[#<0.5,0,1]&/@sData;
tMax=60;
tPlot=Round[tMax/tfinalSim Length[trajectory]]; (** should be at tMax s **)
position=trajectory[[tPlot,1;;3]];
orientation=trajectory[[tPlot,4;;6]];
{front,end}=Transpose[(1+rEcoli)position+0.5Transpose[2rEcoli{#,-#}&@orientation]];
{runColor,tumbleColor}={Darker@Gray,Red};
plotTrajectory=Show[
	Graphics3D[{Opacity[0.1],Sphere[]}],
	Graphics3D[{ColorData[97,1],AbsoluteThickness[3],Line[#&/@trajectory[[;;tPlot,1;;3]],VertexColors->(If[#==0,runColor,tumbleColor]&/@rtData)]}],
	Graphics3D[{If[rtData[[tPlot]]==0,runColor,tumbleColor],CapsuleShape[{front,end},rEcoli]}],
	Graphics3D[{Black,Arrow[Tube[1.3{{0,0,0},#}]]&/@{{1,0,0},{0,1,0},{0,0,1}}}]
,Boxed->False,Axes->False,ImageSize->800,PlotRange->1.4{{-1,1},{-1,1},{-1,1}}
,ViewVertical->{0.376,0.248,0.893},ViewPoint->{2.813,1.728,0.742}
];
If[plotPiecesQ\[And]pthout!="",Export[FileNameJoin[{pthout,"trajectory"<>parameterString<>".png"}],plotTrajectory]];

(** combined plot **)
plot=Grid[{{plotSpectra,Column[{plotData,plotTrajectory}]}}];
If[pthout!="",
	Export[FileNameJoin[{pthout,"combined"<>parameterString<>".png"}],plot]
];

If[debugPlotQ,
	Print[plot];
];

];


(* ::Subsection::Closed:: *)
(*loadDataFromFile*)


(* ::Input:: *)
(*trajectories=loadDataFromFile[pthFileName,nVars,nSteps];*)


loadDataFromFile[pthFileName_,nVars_,nSteps_]:=Module[{type,fstream,rawdata,cdata,trajectories},

fstream=OpenRead[pthFileName,BinaryFormat->True];

(** read data **)
type="Real64";
rawdata=BinaryReadList[fstream,type];
Close[fstream];
cdata=Partition[rawdata,nVars];
trajectories=Partition[cdata,nSteps];

trajectories
]


(* ::Subsection::Closed:: *)
(*readSimulationOutput*)


(* ::Input:: *)
(*{nSteps,dt,times,trajectories,stepsTransient}=readSimulationOutput[tfinalSim,tTransient,stepSize,FileNameJoin[{NotebookDirectory[],"output_figures","runTumbleOutput_young.bin"}],7,10,0.001];*)


(* ::Subsubsection:: *)
(*code*)


Clear[readSimulationOutput];
readSimulationOutput[tfinalSim_,tTransient_,stepSize_:1,pthFnInput_:FileNameJoin[{NotebookDirectory[],"build","runTumbleOutput.bin"}],nVars_:7,nSaveSteps_:10,dtNumeric_:0.001]:=Module[
{pthFileName,nSteps,dt,times,trajectories,dtFine,nStepsFine,timesFine,trajectoriesFine,stepsTransient},
dtFine=dtNumeric nSaveSteps;
nStepsFine=Round[tfinalSim/dtFine]+1;
timesFine=Range[0,tfinalSim,dtFine];
trajectoriesFine=loadDataFromFile[pthFnInput,nVars,nStepsFine];

(** coarsening data in time; Note: Round[nStepsFine/stepSize] might be problematic **)
{nSteps,dt,times,trajectories}={Round[nStepsFine/stepSize],stepSize dtFine,timesFine[[;;;;stepSize]],trajectoriesFine[[All,;;;;stepSize]]};
stepsTransient=Round[tTransient/dt];

{nSteps,dt,times,trajectories,stepsTransient}
]


(* ::Subsection::Closed:: *)
(*runSimulation*)


(* ::Text:: *)
(*TODO: check for errors from execution!*)


(* ::Input:: *)
(*args={u0Dim,diffRotDim,diffS,timeFactor,rDropletDim,asymmetry,swimMultiplier,tfinalSim,nParallelSystems,noiseSeed};*)
(*runSimulation[args];*)


(* ::Subsubsection:: *)
(*code*)


(* ::Input:: *)
(*runSimulation[args]*)


runSimulation[args_]:=Module[{exe,pthExe,arguments,process},
	exe=If[$OperatingSystem=="Windows","run_sim_win.exe","run_sim_unix.exe"];
	pthExe=FileNameJoin[{NotebookDirectory[],"build",exe}];
	If[!FileExistsQ@pthExe,Print["Error: Executable ("<>pthExe<>") does not exist!"]];
	arguments=ToString[DecimalForm[#]]&/@args;
	process=StartProcess[Join[{pthExe},arguments]];
	While[ProcessStatus[process]=="Running",Pause[0.1];];   (** stop async mathematica execution until external process/executable is done **)
]


(* ::Subsection::Closed:: *)
(*runSimulationDebug*)


(* ::Text:: *)
(*TODO: check for errors from execution!*)


(* ::Input:: *)
(*args={u0Dim,diffRotDim,diffS,timeFactor,rDropletDim,asymmetry,swimMultiplier,tfinalSim,nParallelSystems,noiseSeed};*)
(*runSimulationDebug[args];*)


(* ::Subsubsection::Closed:: *)
(*code*)


(* ::Text:: *)
(*shows live window output, Run[] pauses mathematica execution until external program is done*)


(* ::Input:: *)
(*runSimulationDebug[args]*)


runSimulationDebug[args_]:=Module[{exe,pthExe,arguments,cmdString},
	exe=If[$OperatingSystem=="Windows","run_sim_win.exe","run_sim_unix.exe"];
	pthExe=FileNameJoin[{NotebookDirectory[],"build",exe}];
	If[!FileExistsQ@pthExe,Print["Error: Executable ("<>pthExe<>") does not exist!"]];
	arguments=StringJoin@@Riffle[ToString[DecimalForm[#]]&/@args," "];
	cmdString="\""<>pthExe<>"\""<>" "<>arguments; (** for windows cmd **)
	Run[cmdString];
]


(* ::Subsection::Closed:: *)
(*thetaNoiseReconstruction*)


thetaNoiseReconstruction[thetaData_,dt_,gravFactor_]:=Module[{thetaDotData,thetaGravData},
thetaDotData=Differences[thetaData]/dt;
thetaGravData=gravFactor Sin[thetaData[[;;-2]]];

thetaDotData-thetaGravData  (** theta noise **)
]


(* ::Subsection:: *)
(*updateIntensityLookupMap*)


(* ::Text:: *)
(*updates global variable thetaPhiIntensityInterpolation*)


(* ::Input:: *)
(*path=FileNameJoin[{NotebookDirectory[],"..","Droplet optics","prop_test_dropletTilt_thetaPhi_analyticDirector_wvl_470_DnLC_0.20","zDistance_0.0"}];*)
(*path=FileNameJoin[{NotebookDirectory[],"..","Droplet optics","RGB"}]*)
(*thetaPhiIntensityInterpolation=updateIntensityLookupMap[path];*)


(* ::Code:: *)
(*rawString="0,630226339\t0,650020406\t0,689197473\t0,717853051\t0,713290728\t0,654277047\t0,566341708\t0,455507696\t0,309285571\t0,128632962\t0,239113181\t0,733882383\t0,339818443\t0,034834678\t0\t0,068419135*)
(*0,630226339\t0,650152628\t0,689355223\t0,719570958\t0,717274742\t0,660631244\t0,577113895\t0,468189579\t0,316722419\t0,129818381\t0,247327654\t0,751471873\t0,363761802\t0,058072415\t0,016610914\t0,080852938*)
(*0,630226339\t0,650029243\t0,688500033\t0,72233748\t0,727031044\t0,680592552\t0,608075181\t0,497458485\t0,33418296\t0,135175019\t0,260609119\t0,803221183\t0,438199669\t0,119015428\t0,07870989\t0,131495051*)
(*0,630226339\t0,649701305\t0,687974744\t0,726963296\t0,743272128\t0,714940561\t0,659793072\t0,553020313\t0,362899412\t0,147901741\t0,277911909\t0,885610387\t0,555869311\t0,218768952\t0,157403089\t0,213451587*)
(*0,630226339\t0,64955501\t0,687408545\t0,732298334\t0,762349446\t0,75555309\t0,722216386\t0,622943641\t0,402085643\t0,165674509\t0,295892175\t0,96331061\t0,697705845\t0,341970655\t0,242267861\t0,2940499*)
(*0,630226339\t0,649302675\t0,686636485\t0,738172407\t0,781243158\t0,797500279\t0,787409821\t0,698130004\t0,450949628\t0,18348426\t0,303611468\t0,953760169\t0,83785851\t0,452647571\t0,321843741\t0,386442831*)
(*0,630226339\t0,649055576\t0,685370555\t0,742570434\t0,797122922\t0,833614633\t0,840085342\t0,767275707\t0,501594358\t0,19511982\t0,215505614\t0,827300038\t0,951890075\t0,537451468\t0,369737987\t0,423047455*)
(*0,630226339\t0,648694256\t0,684327831\t0,742554724\t0,802617346\t0,849894009\t0,867675614\t0,816328534\t0,543798137\t0,186121178\t0,066531368\t0,629502717\t1\t0,577370485\t0,390194475\t0,445857106";*)
(**)
(*xlsxData=Partition[ToExpression[StringJoin[{"{",StringReplace[rawString,{","->".","\t"->",","\n"->","}],"}"}]],16]\[Transpose];*)
(*Export[FileNameJoin[{$HomeDirectory,"Desktop","normIntegratedIntensity.mat"}],xlsxData]*)


(* ::Code:: *)
(*data1=Import[FileNameJoin[{$HomeDirectory,"Desktop","normIntegratedIntensity.mat"}]][[1]]//Dimensions*)
(*data2=Import[FileNameJoin[{$HomeDirectory,"Desktop","pathogen_sensing","Data","numerical_droplet_optics","lookup_map","normIntegratedIntensity.mat"}]][[1]]//Dimensions*)


(* ::Subsubsection:: *)
(*create normIntegratedIntensity.mat from data*)


(* ::Subsubsection::Closed:: *)
(*fit low-order trig function*)


(* ::Input:: *)
(*Plot[Sin[2  \[Pi]/2 (\[Theta]/(\[Pi]/2))^0.4]^2(1-0.3Cos[2 0]^2),{\[Theta],0,\[Pi]/2},ImageSize->300,Evaluate[ops],FrameLabel->{"polar angle \!\(\**)
(*StyleBox[\"\[Theta]\",\nFontSlant->\"Italic\"]\)","intensity transfer function \!\(\**)
(*StyleBox[\"I\",\nFontSlant->\"Italic\"]\)\!\(\**)
(*StyleBox[\"(\",\nFontSlant->\"Italic\"]\)\!\(\**)
(*StyleBox[\"\[Theta]\",\nFontSlant->\"Italic\"]\)\!\(\**)
(*StyleBox[\")\",\nFontSlant->\"Italic\"]\)"}]*)


(* ::Input:: *)
(*ops={Frame->True,FrameStyle->Directive[Black,18,AbsoluteThickness[2]],PlotRange->All,AspectRatio->1,ImageSize->{Automatic,400},ImagePadding->{{60,10},{60,5}}};*)
(*f[t_]:=(*0.5*)\[Pi]/2/2.1(Sin[t]+1.1)*)
(*p1=Plot[f[x],{x,0,5 2\[Pi]},Evaluate[ops],FrameLabel->{"t","\[Theta](t)"}];*)
(*p2=Plot[fit\[Theta]\[Phi]1[f[x],0],{x,0,5 2\[Pi]},Evaluate[ops],FrameLabel->{"t","I(\[Theta](t))"}];*)
(*fData1=Table[f[x],{x,Subdivide[0,5. 2\[Pi],1000]}];*)
(*fData2=Table[fit\[Theta]\[Phi]1[f[x],0],{x,Subdivide[0,5. 2\[Pi],1000]}];*)
(*p3=ListLinePlot[Abs[Fourier[fData1]][[2;;50]]^2,Evaluate[ops],FrameLabel->{"frequency index \!\(\*SubscriptBox[\(f\), \(i\)]\)","PSD[I(\[Theta](t))]"}];*)
(*p4=ListLinePlot[Abs[Fourier[fData2]][[2;;50]]^2,Evaluate[ops],FrameLabel->{"frequency index \!\(\*SubscriptBox[\(f\), \(i\)]\)","PSD[I(\[Theta](t))]"}];*)
(**)
(*plot=Grid[{{p1,p3},{p2,p4}}]*)
(*Export[FileNameJoin[{$HomeDirectory,"Desktop","Ift_test_osci_large.png"}],plot]*)


(* ::Input:: *)
(*p1=Plot3D[thetaPhiIntensityInterpolation[\[Pi]/180 \[Theta],\[Pi]/180 \[Phi]],{\[Theta],0,90},{\[Phi],0,45},BoxRatios->{1,0.5,0.5}];*)
(*(*\[Pi]/180 \[Phi]]*)*)
(*fit\[Theta]\[Phi]1[\[Theta]_,\[Phi]_]:=Sin[2  \[Pi]/2 (\[Theta]/(\[Pi]/2))^0.4]^2(1-0.3Cos[2\[Phi]]^2)*)
(*fit\[Theta]\[Phi]2[\[Theta]_,\[Phi]_]:=Sin[\[Pi] (\[Theta]/(\[Pi]/2))^0.3]^2 (1-0.3Cos[2\[Phi]]^2)*)
(*fit\[Theta]\[Phi]3[\[Theta]_,\[Phi]_]:=\[Theta] Exp[-10\[Theta]] (1-0.3Cos[2\[Phi]]^2)*)
(*p2=Plot3D[fit\[Theta]\[Phi]1[\[Pi]/180 \[Theta],\[Pi]/180\[Phi]],{\[Theta],0,90},{\[Phi],0,45},BoxRatios->{1,0.5,0.5},PlotRange->All];*)
(*p3=Plot3D[fit\[Theta]\[Phi]2[\[Pi]/180 \[Theta],\[Pi]/180\[Phi]],{\[Theta],0,90},{\[Phi],0,45},BoxRatios->{1,0.5,0.5},PlotRange->All,PlotPoints->20];*)
(*p4=Plot3D[fit\[Theta]\[Phi]3[\[Pi]/180 \[Theta],\[Pi]/180\[Phi]],{\[Theta],0,90},{\[Phi],0,45},BoxRatios->{1,0.5,0.5},PlotRange->All,PlotPoints->20];*)
(*{p1,p2,p3,p4}*)


(* ::Subsubsection::Closed:: *)
(*intensity map visualization*)


(* ::Input:: *)
(*name=FileNameJoin[{"prop_test_dropletTilt_thetaPhi_analyticDirector_wvl_470_DnLC_0.20","zDistance_0.0"}];*)
(*nameLookupMap="prop_true_compiled_test_dropletTilt_thetaPhi_analyticDirector_wvl_550_DnLC_0.18\\zDistance_0.0"*)
(*updateIntensityLookupMap[name];*)


(* ::Input:: *)
(*Plot3D[thetaPhiIntensityInterpolation[\[Pi]/180 \[Theta],\[Pi]/180 \[Phi]],{\[Theta],0,90},{\[Phi],0,45},BoxRatios->{1,0.5,0.5}]*)
(**)
(*ContourPlot[thetaPhiIntensityInterpolation[\[Pi]/180 \[Theta],\[Pi]/180 \[Phi]],{\[Theta],0,90},{\[Phi],0,45},AspectRatio->1/2,PlotRangePadding->None,FrameLabel->{"\[Theta] (\[Degree])","\[Phi] (\[Degree])"},FrameStyle->Directive[Black,16]]*)


(* ::Input:: *)
(*(** \[Theta]=radial coordinate, \[Phi]=angle coordinate **)*)
(*data=Flatten[Table[{\[Theta] Cos[\[Phi]],\[Theta] Sin[\[Phi]],thetaPhiIntensityInterpolation[-Abs[\[Theta]-\[Pi]/2.0]+\[Pi]/2.0,-Abs[Mod[\[Phi],\[Pi]/2.0]-\[Pi]/4.0]+\[Pi]/4.0]},{\[Theta],Subdivide[0,\[Pi]/2,100-1]},{\[Phi],Most[Subdivide[0,2\[Pi],100-1]]}],1];*)
(*ListDensityPlot[data,PlotRange->All,ColorFunction->GrayLevel,Epilog->{Darker@Red,*)
(*Arrow[{{0,0},{\[Pi]/2,0}}],Text[Style["\[Theta]",18],Scaled[{0.92,0.5}]],Arrow@ResourceFunction["SplineCircle"][{0,0},1.1\[Pi]/2,{1,0},{5Degree,90 Degree}],Text[Style["\[Phi]",18],Scaled[{0.83,0.83}]]},Frame->False,PlotRangePadding->Scaled[0.1]]*)
(**)
(*ListPlot3D[data,PlotRange->All,ColorFunction->GrayLevel,Boxed->False,Axes->False]*)


(* ::Subsubsection:: *)
(*code*)


updateIntensityLookupMap[path_]:=Module[{thetaPhiIntensities,\[Theta]s,\[Phi]s,\[Theta]\[Phi]IData},

	(** get externally computed dataset of tilt angles -> intensity **)
	thetaPhiIntensities=Import[FileNameJoin[{path,"normIntegratedIntensity.mat"}]][[1]];
	
	(** generate continuous lookup map form tilt angles to intensity**)
	\[Theta]s=N@Subdivide[0,\[Pi]/2,Dimensions[thetaPhiIntensities][[1]]-1];
	\[Phi]s=N@Subdivide[0,\[Pi]/4,Dimensions[thetaPhiIntensities][[2]]-1];
	\[Theta]\[Phi]IData=Flatten[Table[{\[Theta]s[[i]],\[Phi]s[[j]],thetaPhiIntensities[[i,j]]},{i,Length[\[Theta]s]},{j,Length[\[Phi]s]}],1];
	
	(** output **)
	Interpolation[\[Theta]\[Phi]IData,InterpolationOrder->{1,1}]  
]
