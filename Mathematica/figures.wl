(* ::Package:: *)

(* ::Title:: *)
(*Paper Figures*)


(* ::Chapter::Closed:: *)
(*load custom functions from external package*)


Get[FileNameJoin[{NotebookDirectory[],"functions.wl"}]];
Get[FileNameJoin[{NotebookDirectory[],"droplet_optics.wl"}]];


(* ::Chapter::Closed:: *)
(*Constants*)


(* ::Subsection::Closed:: *)
(*parametersExperiment: parameters*)


parametersExperiment=<|
	"nDays"->6
	,"fps"->30
	,"dt":>1.0/parametersExperiment["fps"]
	,"movieLengths"->{8344, 7252,8046,9225,5244,2251}
	,"numberDroplets"->{71,50,83,124,62,14}
|>;


(* ::Subsection::Closed:: *)
(*parametersFourier: for spectral analysis*)


parametersFourier=<|
	"windowingQ"->True
	,"fMax"->1.0
	,"padFactor"->2
	,"nPaddedSamples":>parametersFourier["padFactor"] Max[parametersExperiment["movieLengths"]]
|>;


(* ::Subsection::Closed:: *)
(*parametersPlot: journal rules for plotting*)


parametersPlot=<|
	"pts" -> 99.6923076923077`  (** weird correction factor for mma 13.3 **)
	,"cm" :> parametersPlot["pts"]/2.54 (* centimetre *)
	,"mm" :> 0.1 parametersPlot["pts"]/2.54 (* millimetre *)
	,"fontSize" :> 10 parametersPlot["pts"]/72
	,"fontFamily" -> "Arial"
	,"lineThickness" :> 1 parametersPlot["pts"]/72
	,"timeLabel"->Row[{"Time ",Style["t",FontSlant->"Italic"]," (s)"}]
	,"intensityLabel"->Row[{"Intensity ",Style["I",FontFamily->"Century",FontSlant->Italic]}]
	,"frequencyLabel"->Row[{"Frequency ",Style["f",FontSlant->"Italic"]," (Hz)"}]
	,"spectrumLabel"->Text[Style["Intensity\nspectrum", LineSpacing -> {0.85, 0}]]
	,"widthLabel"->"Width \!\(\*StyleBox[SubscriptBox[\"f\", \"HM\"],\nFontSlant->\"Italic\"]\)\!\(\*StyleBox[\" \",\nFontSlant->\"Italic\"]\)(Hz)"
	,"ops" :>{
		Axes->None
		,Frame->True
		,FrameStyle->Directive[Black,parametersPlot["fontSize"],AbsoluteThickness[parametersPlot["lineThickness"]],FontFamily:>parametersPlot["fontFamily"]]
	}
|>;


(* ::Subsection::Closed:: *)
(*parametersColor: colors*)


(*{colorYoung,colorOld}={
TemplateBox[Association["color" -> RGBColor[0.5, 0.81, 0]], "RGBColorSwatchTemplate"],};*)
(*{colorYoung,colorOld}={
TemplateBox[Association["color" -> RGBColor[0.5, 0.81, 0]], "RGBColorSwatchTemplate"],};*)
(*{colorYoung,colorOld}={
TemplateBox[Association["color" -> RGBColor[0.5, 0.81, 0]], "RGBColorSwatchTemplate"],
TemplateBox[Association["color" -> RGBColor[0.62, 0.62, 0.62]], "RGBColorSwatchTemplate"]};*)

parametersColor=<|
	"young"->RGBColor[0.5, 0.81, 0]
	,"old"->Lighter@RGBColor[0.3, 0.17, 0.82]
	,"ages":>Table[Blend[{parametersColor["young"],parametersColor["old"]},x],{x,Subdivide[0,1,parametersExperiment["nDays"]-1]}]
	,"youngOld":>(parametersColor[#]&/@{"young","old"})
|>;


(* ::Subsection::Closed:: *)
(*parametersPaths: paths*)


basePth=NotebookDirectory[];
pthData=FileNameJoin[{ParentDirectory[basePth],"Data"}];
pthVideo=FileNameJoin[{ParentDirectory[basePth],"Videos"}];
debugQ=False;

parametersPaths=<|
	"rawMoviesExperiment"->pthVideo
	,"intensityDataExperiment"->FileNameJoin[{pthData,"experimental_intensity_time_traces"}]
	,"intensitySpectraExperiment"->FileNameJoin[{pthData,"experimental_intensity_spectra"}]
	,"dropletTrajectoriesExperiment"->FileNameJoin[{pthData,"experimental_droplet_trajectories"}]
	,"widthsExperiment"->FileNameJoin[{pthData,"experimental_widths"}]
	,"dropletOpticsNumerical"->FileNameJoin[{pthData,"numerical_droplet_optics"}]
	,"runTumbleTrajectoriesNumerical"->FileNameJoin[{pthData,"numerical_runTumble_trajectories"}]
	,"intensitySpectraNumerical"->FileNameJoin[{pthData,"numerical_intensity_spectra"}]
	,"parameterSweepNumerical"->FileNameJoin[{pthData,"numerical_parameter_sweep"}]
	,"figs"->If[debugQ,
		FileNameJoin[{basePth,"Figures","Figures debug"}],
		FileNameJoin[{basePth,"Figures","Figures current"}]
	]
|>;


(* ::Subsection::Closed:: *)
(*p: all parameters together*)


(* ::Text:: *)
(*Example call:*)


(* ::Input:: *)
(*p["experiment","dt"]*)
(*p["color","youngOld"]*)
(*p["path","figs"]*)
(*p@@{"path","rawMoviesExperiment"}*)


p=<|
	"experiment"->parametersExperiment
	,"fourier"->parametersFourier
	,"plot"->parametersPlot
	,"color"->parametersColor
	,"path"->parametersPaths
|>;


(* ::Chapter:: *)
(*General functions*)


(* ::Subsection::Closed:: *)
(*colormaps: ageColormap, timeColormap*)


ageColormap[x_]:=Blend[p["color",#]&/@{"young","old"},x];   (** bacterium age **)
timeColormap[x_]:=Blend[{RGBColor[0.9294117647058824, 0.12156862745098039`, 0.1411764705882353],RGBColor[1., 0.7843137254901961, 0.0392156862745098]},x];                          (** experiment time **)


(* ::Subsection::Closed:: *)
(*calculateDailyMeanSpectra*)


(* ::Input:: *)
(*{allIntensities,freqs,meanSpectra,stdSpectra,allSpectra,weightedMeanSpectra,weightedStdSpectra,allWeightedSpectra}=calculateDailyMeanSpectra[p];*)


(* ::Text:: *)
(*weighted spectra account for lengths of timeseries with factor dt/L (timestep/duration). This is important for correct weighting of spectra when calculating the mean spectrum. Derivation for prefactor is in another notebook.*)


(* ::Subsubsection:: *)
(*test: show single spectra from experiments*)


(* ::Input:: *)
(*individualSpectra=allSpectra[[3]];*)
(*Manipulate[*)
(*	ListLinePlot[{{freqs,individualSpectra[[i]]}\[Transpose],{freqs,meanSpectraEx[[3]]}\[Transpose]},PlotRange->{{0,1},All},Frame->True,AspectRatio->1]*)
(*,{i,1,Length[individualSpectra],1,Appearance->"Open"}]*)


(* ::Subsubsection::Closed:: *)
(*code*)


calculateDailyMeanSpectra[p_]:=Module[
{dayIndices,allIntensities,meanSpectra,stdSpectra,allSpectra,intensityParts,times,paddedMemory,powerSpectralDensities,
freqs,spectra,meanSpectrum,stdSpectrum,weightedMeanSpectra,allWeightedSpectra,lengths,weightedSpectra,weightedMeanSpectrum,weightedStdSpectrum,weightedStdSpectra},
	
	(** calculate mean spectrum for each day **)
	dayIndices=Range[p["experiment","nDays"]];
	Print["Runtime - calculateDailyMeanSpectra: ",AbsoluteTiming[
		{allIntensities,meanSpectra,stdSpectra,allSpectra,weightedMeanSpectra,weightedStdSpectra,allWeightedSpectra}=Table[
		
			(** get intensity **)
			intensityParts=getIntensitySingleFromExperiment[p,dayIndex];
			lengths=Length/@intensityParts;
			
			(** get mean/std spectrum **)
			times=Range[0,Length[intensityParts[[1]]]-1]p["experiment","dt"];
			paddedMemory=ConstantArray[0.0,p["fourier","nPaddedSamples"]];
			powerSpectralDensities=calculateTemporalFreqSpectrum[p["experiment","dt"],#,p["fourier","windowingQ"],p["fourier","nPaddedSamples"],paddedMemory][[{1,2}]]\[Transpose]&/@intensityParts;
			freqs=powerSpectralDensities[[1,All,1]];
			spectra=powerSpectralDensities[[All,All,2]];
			weightedSpectra=MapThread[#1 #2&,{p["experiment","dt"]/lengths,spectra}];
			meanSpectrum=Mean[spectra];
			stdSpectrum=StandardDeviation[spectra];
			weightedMeanSpectrum=Mean[weightedSpectra];
			weightedStdSpectrum=StandardDeviation[weightedSpectra];

			(** output **)
			{intensityParts,meanSpectrum,stdSpectrum,spectra,weightedMeanSpectrum,weightedStdSpectrum,weightedSpectra}
		,{dayIndex,dayIndices}]\[Transpose];
	][[1]],"\[ThinSpace]s"];
	
	(** output **)
	{allIntensities,freqs,meanSpectra,stdSpectra,allSpectra,weightedMeanSpectra,weightedStdSpectra,allWeightedSpectra}
]


(* ::Subsection::Closed:: *)
(*filterBadDroplets*)


(* ::Input:: *)
(*goodDropletIndices = filterBadDroplets[p["numberDroplets"]];*)


(* ::Subsubsection::Closed:: *)
(*code*)


filterBadDroplets[numberDroplets_]:=Module[{badLists},

(** filter out bad droplets due to irregularities found by eye **)
badLists={
	{1,3,7,9,10,16,20,23,27,29,41,43,44,46,49,50,54,55,59,64},
	{4,5,17,18,21,22,23,33,36,38,47,49},
	{4,6,7,9,15,16,21,24,36,39,40,42,43,45,59,62,66,71,73},
	{3,6,7,11,13,14,15,24,30,31,32,33,40,43,44,51,57,65,77,79,80,82,85,92,96,98,103,105,108,114,117,118,124},
	{5,13,21,32,33,34,52,55},
	{}
};

(** output: lists of good droplets **)
MapThread[Complement[Range[#2],#1]&,{badLists,numberDroplets}]
]


(* ::Subsection::Closed:: *)
(*getDropletFramesFromMovie*)


(* ::Input:: *)
(*singleDropletImages=getDropletFramesFromMovie[nSingleDropletImages,p,{timeStart,timeDuration,dayIndex,dropletIndices[[dayIndex]]}]*)


(* ::Subsubsection:: *)
(*debug: test images vs trajectory with Manipulates*)


(* ::Input:: *)
(*dayIndex=1;*)
(*dropletIndex=2;*)
(*{timeStart,timeDuration}={10,70};*)
(*nSingleDropletImages=4;*)
(**)
(*(** clean and sort droplet trajectories in the same order as imported intensities **)*)
(*trajectories0=getTrajectories[p,dayIndex];*)
(*{trajectories,starts}=sortCleanTrajectories[trajectories0];*)
(**)
(*(** get specific droplet trajectory **)*)
(*trajectory=trajectories[[dropletIndex]];*)
(**)
(*(** get frame indices for sparse sampling of video **)*)
(*trajectoryStart=starts[[dropletIndex]];*)
(*{frameStart,frameEnd}=trajectoryStart+Round[{timeStart,timeDuration}/p["experiment","dt"]];*)
(*frameIndices=Round[Subdivide[frameStart,frameEnd,nSingleDropletImages-1]];  *)
(**)
(*(** get video frames **)*)
(*movieFile=FileNameJoin[{p["path","rawMoviesExperiment"],"Day "<>ToString[dayIndex]<>".avi"}];*)
(*(*frameIndices=Round[Subdivide[frameStart,frameEnd,nSingleDropletImages-1]];*)*)
(*frameIndices=Range[frameStart,frameEnd];*)
(*{{nx,ny},dropletImages}=getMovieFrames[movieFile,frameIndices];*)


(* ::Input:: *)
(*dropletImages//Length*)


(* ::Input:: *)
(*dropletImages[[1]]//ImageDimensions*)


(* ::Input:: *)
(*dropletImages[[1]]*)


(* ::Input:: *)
(*Graphics[{Rectangle[{0,0},{720,480}]}]*)


(* ::Input:: *)
(*boxLength=16;*)
(*halfBoxLength=Round[0.5boxLength];*)
(*Manipulate[*)
(*	position=trajectory[[frameIndices[[i]]-trajectoryStart+1]];*)
(*	boxCorners=Round[{{Max[#[[1]]-halfBoxLength,1],Min[#[[1]]+halfBoxLength,nx]},{Max[#[[2]]-halfBoxLength,1],Min[#[[2]]+halfBoxLength,ny]}}]&@position;*)
(*	img=ImageTake[dropletImages[[i]],ny-#[[2]],#[[1]]]&@boxCorners;*)
(*Show[img,PlotLabel->position]*)
(*,{i,1,Length[frameIndices],1}]*)
(**)


(* ::Input:: *)
(*boxLength=16;*)
(*halfBoxLength=Round[0.5boxLength];*)
(*Manipulate[*)
(*	p1=Show[*)
(*		dropletImages[[frame-frameStart+1]],*)
(*		ListPlot[trajectories[[dropletIndex]]]*)
(*	];*)
(*	position=trajectory[[frame-trajectoryStart+1]];*)
(*	boxCorners=Round[{{Max[#[[1]]-halfBoxLength,1],Min[#[[1]]+halfBoxLength,nx]},{Max[#[[2]]-halfBoxLength,1],Min[#[[2]]+halfBoxLength,ny]}}]&@position;*)
(*p2=ImageTake[dropletImages[[frame-frameStart+1]],ny-#[[2]],#[[1]]]&@boxCorners;*)
(*Grid[{{p1,p2}}]*)
(*,{{dropletIndex,16},1,Length[trajectories],1,Appearance->"Open"}*)
(*,{frame,frameStart,frameEnd,1,Appearance->"Open"}*)
(*,TrackedSymbols:>{dropletIndex,frame}]*)


(* ::Subsubsection:: *)
(*code*)


getDropletFramesFromMovie[nSingleDropletImages_,p_,{timeStart_,timeDuration_,dayIndex_,dropletIndex_}]:=Module[{trajectoryStart,frameStart,frameEnd,frameIndices,
movieFile,nx,ny,dropletImages,trajectories0,trajectories,starts,trajectory,start,boxLength,halfBoxLength,singleDropletImages,position,boxCorners,filterLists,dropletIndexGood},
	
	(** clean and sort droplet trajectories in the same order as imported intensities **)
	trajectories0=getTrajectories[p,dayIndex];
	{trajectories,starts}=sortCleanTrajectories[trajectories0];
	
	(** filter out bad droplets from manual observation **)
	filterLists = filterBadDroplets[p["experiment","numberDroplets"]];
	dropletIndexGood=filterLists[[dayIndex,dropletIndex]];
	
	(** get specific droplet trajectory **)
	trajectory=trajectories[[dropletIndexGood]];
	
	(** get frame indices for sparse sampling of video **)
	trajectoryStart=starts[[dropletIndexGood]];
	{frameStart,frameEnd}=trajectoryStart+Round[{timeStart,timeDuration}/p["experiment","dt"]];
	frameIndices=Round[Subdivide[frameStart,frameEnd,nSingleDropletImages-1]];  
	
	(** get video frames **)
	movieFile=FileNameJoin[{p["path","rawMoviesExperiment"],"Day "<>ToString[dayIndex]<>".avi"}];
	{{nx,ny},dropletImages}=getMovieFrames[movieFile,frameIndices];
	
	(** DEBUG **)
	Print[
		"frame indices: ",frameIndices,"\n"
		,"positions (x,y): ",trajectory[[#]]&/@(frameIndices-trajectoryStart+1)
	];
	
	(** get droplet frames **)
	boxLength=16;
	halfBoxLength=Round[0.5boxLength];
	singleDropletImages=Table[
		position=trajectory[[frameIndices[[i]]-trajectoryStart+1]];
		boxCorners=Round[{{Max[#[[1]]-halfBoxLength,1],Min[#[[1]]+halfBoxLength,nx]},{Max[#[[2]]-halfBoxLength,1],Min[#[[2]]+halfBoxLength,ny]}}]&@position;
		ImageTake[dropletImages[[i]],ny-#[[2]],#[[1]]]&@boxCorners
	,{i,Length[frameIndices]}];
	
	
	(** output **)
	singleDropletImages
]


(* ::Subsection::Closed:: *)
(*getIntensityImages*)


(* ::Input:: *)
(*getIntensityImages[thetaData,phiData,{stepStart,stepEnd},p,nIntensityImages]*)


(* ::Text:: *)
(*Image gets rotated appropriately if outside of \[Phi] symmetry bounds!*)


(* ::Subsubsection::Closed:: *)
(*debug: read out image variables*)


(* ::Input:: *)
(*pthIntensityImages="C:\\Users\\Jan\\Dropbox (MIT)\\Bacteria Droplet Manuscript\\Code\\Droplet optics\\rgb avg for fig3\\RGB_averaged_images";*)
(*phiThetaValues=ToExpression[StringSplit[FileBaseName[#],"_"][[{4,6}]]]&/@FileNames[All,pthIntensityImages];*)
(*Tally[phiThetaValues[[All,1]]]*)
(*Tally[phiThetaValues[[All,2]]]*)


(* ::Subsubsection::Closed:: *)
(*debug: image cropping*)


(* ::Input:: *)
(*Import[FileNameJoin[{pthIntensityImages,"afterDropletIntensityXY_RGB_theta_0.000_phi_0.000.png"}]]*)


(* ::Subsubsection::Closed:: *)
(*debug: intensity*)


(* ::Input:: *)
(*MinMax[Table[thetaPhiIntensityInterpolation[\[Theta],\[Phi]],{\[Theta],Subdivide[0,\[Pi]/2,100-1]},{\[Phi],Subdivide[0,\[Pi]/4,100-1]}]]*)


(* ::Input:: *)
(*plot=DensityPlot[thetaPhiIntensityInterpolation[\[Theta],\[Phi]],{\[Theta],0,\[Pi]/2},{\[Phi],0,\[Pi]/4},PlotLegends->Automatic,FrameLabel->{"\[Theta]","\[Phi]"},PlotRangePadding->None,FrameStyle->Directive[Black,AbsoluteThickness[lineThickness],fontSize],ColorFunction->"BlueGreenYellow",AspectRatio->1/2]*)


(* ::Subsubsection::Closed:: *)
(*debug: rotate image based on phi*)


(* ::Input:: *)
(*phiData=Table[x,{x,Subdivide[0,2 2\[Pi],100-1]}];*)
(*thetaData=ConstantArray[0.49\[Pi],100];*)
(*{stepStart,stepEnd}={1,100};*)
(*nIntensityImages=100;*)
(*images=getIntensityImages[thetaData,phiData,{stepStart,stepEnd},p,nIntensityImages];*)


(* ::Input:: *)
(*Manipulate[*)
(*n=Quotient[phiData[[i]],\[Pi]/4];*)
(*Grid[{{{N@phiData[[i]],n},SpanFromLeft},{ListPlot[Quotient[phiData,\[Pi]/4],Epilog->{Red,InfiniteLine[{{#,0},{#,1}}&@i]},ImageSize->{Automatic,300},Frame->True,AspectRatio->1,Axes->False,PlotRange->All],images[[i]],If[EvenQ[n],ImageRotate[#,-\[Pi]/4 n],ImageRotate[ImageReflect[#,Right->Top],-\[Pi]/4 (n-1)]]&@images[[i]]}}]*)
(*,{i,1,nIntensityImages,1,Appearance->"Open"}]*)


(* ::Subsubsection:: *)
(*code*)


Clear[getIntensityImages];
getIntensityImages[thetaData_,phiData_,{stepStart_,stepEnd_},p_,nImages_:5]:=Module[
{n\[Theta],n\[Phi],nAllImages,timeSteps,thetaPhiTrajPoints,thetaPhiSamplePoints,thetaPhiTrajStrings,images,fileNames,desiredFileNames,pthNumericalIntensityImages,nSymmetryJumps,prefix},

	pthNumericalIntensityImages=FileNameJoin[{p["path","dropletOpticsNumerical"],"intensity_images"}];
	Print["Source for intensity images: ",pthNumericalIntensityImages];
	fileNames=FileNames[All,pthNumericalIntensityImages];
	prefix=StringSplit[FileBaseName[fileNames[[1]]],"_theta_"][[1]];
	{n\[Theta],n\[Phi]}=Length/@DeleteDuplicates/@Transpose[StringSplit[FileBaseName[#],"_"][[{-3,-1}]]&/@fileNames];
	nAllImages=Length[FileNames["*.png",pthNumericalIntensityImages]];
	If[n\[Theta] n\[Phi]!=nAllImages,Print[Style["Warning: Less than expected number of images ("<>ToString[n\[Theta] n\[Phi]]<>") were found in the given image directory ("<>ToString[nAllImages]<>")!",Orange]];];
	
	(** get files in useful order **)
	timeSteps=Round[Subdivide[stepStart,stepEnd,nImages-1]];
	thetaPhiTrajPoints={-Abs[#-\[Pi]/2.0]+\[Pi]/2.0&/@thetaData[[timeSteps]],-Abs[Mod[#,\[Pi]/2.0]-\[Pi]/4.0]+\[Pi]/4.0&/@phiData[[timeSteps]]}\[Transpose];
	thetaPhiSamplePoints=N@Flatten[Table[{\[Theta],\[Phi]},{\[Theta],Subdivide[0,\[Pi]/2,n\[Theta]-1]},{\[Phi],Subdivide[0,\[Pi]/4,n\[Phi]-1]}],1];
	thetaPhiTrajStrings=ToString[DecimalForm[#,{4,3}]]&/@Nearest[thetaPhiSamplePoints,#][[1]]&/@thetaPhiTrajPoints;
	desiredFileNames=prefix<>"_theta_"<>#[[1]]<>"_phi_"<>#[[2]]<>".png"&/@thetaPhiTrajStrings;
	images=Import[FileNameJoin[{pthNumericalIntensityImages,#}]]&/@desiredFileNames;
	
	(** output: rotated images **)
	nSymmetryJumps=Quotient[phiData[[timeSteps]],\[Pi]/4];
	MapThread[If[EvenQ[#1],ImageRotate[#2,-\[Pi]/4 #1],ImageRotate[ImageReflect[#2,Right->Top],-\[Pi]/4 (#1-1)]]&,{nSymmetryJumps,images}]
]


(* ::Subsection::Closed:: *)
(*getIntensitiesFromExperiment*)


(* ::Text:: *)
(*get intensities of all droplets for all days*)


(* ::Code:: *)
(*dayIndices={1,6};*)
(*allIntensities=getIntensitiesFromExperiment[p,dayIndices];*)
(*allIntensities[[1]]//Dimensions*)


(* ::Subsubsection:: *)
(*code*)


Clear[getIntensitiesFromExperiment];
getIntensitiesFromExperiment[p_,dayIndices_]:=Module[{filterQ,filterLists,intensityParts},
	filterQ=True;
	
	Table[
		intensityParts=Import[FileNameJoin[{p["path","intensityDataExperiment"],"all_intensityParts_"<>ToString[dayIndex]<>".dat"}]];
		(** output **)
		intensityParts = If[filterQ,
			filterLists = filterBadDroplets[p["experiment","numberDroplets"]];
			intensityParts[[filterLists[[dayIndex]]]]
		,
			intensityParts
		]
	,{dayIndex,dayIndices}]
]


(* ::Subsection::Closed:: *)
(*getIntensitySingleFromExperiment*)


(* ::Text:: *)
(*get intensities of all droplets for a single day*)


(* ::Code:: *)
(*dayIndex=1;*)
(*intensity=getIntensitySingleFromExperiment[p,dayIndex];*)


(* ::Subsubsection:: *)
(*code*)


Clear[getIntensitySingleFromExperiment];
getIntensitySingleFromExperiment[p_,dayIndex_]:=Module[{filterQ,filterLists,intensityParts},
	filterQ=True;
	
	intensityParts=Import[FileNameJoin[{p["path","intensityDataExperiment"],"all_intensityParts_"<>ToString[dayIndex]<>".dat"}]];
	intensityParts=If[filterQ,
		filterLists = filterBadDroplets[p["experiment","numberDroplets"]];
		intensityParts[[filterLists[[dayIndex]]]]
	,
		intensityParts
	];
	
	intensityParts
]


(* ::Subsection::Closed:: *)
(*getMovieFrames*)


(* ::Subsubsection:: *)
(*code*)


getMovieFrames[movieFile_,frameIndices_]:=Module[{nx,ny,dropletImages},
	
	(** get info about movie **)
	{nx,ny}=Import[movieFile,{"AVI","ImageSize"}];
	
	(** get certain frames **)
	dropletImages=ColorConvert[#,"Grayscale"]&/@Import[movieFile,{"AVI","ImageList",frameIndices}];
	
	(** output **)
	{{nx,ny},dropletImages}
]


(* ::Subsection::Closed:: *)
(*getTrajectories*)


getTrajectories[p_,dayIndex_]:=Import[FileNameJoin[{p["path","dropletTrajectoriesExperiment"],"trajectories_"<>ToString[dayIndex]<>".hdf5"}],{"Datasets","Dataset1"}]


(* ::Subsection::Closed:: *)
(*getWidthAngleMapping*)


(* ::Input:: *)
(*getWidthAngleMapping[parameterPointsStartEnd,widthInterpolationLogLogLin,angleInterpolationLogLogLin]*)


getWidthAngleMapping[parameterPointsStartEnd_,widthInterpolationLogLogLin_,angleInterpolationLogLogLin_]:=Module[{intersectionPointLogLogLinear,width,meanAngleInterpolation},
Table[
	intersectionPointLogLogLinear=lineLinear[t,parameterPointsStartEnd[[1]],parameterPointsStartEnd[[2]]];
	width=widthInterpolationLogLogLin@@intersectionPointLogLogLinear;
	meanAngleInterpolation=angleInterpolationLogLogLin@@intersectionPointLogLogLinear;
	{width,meanAngleInterpolation}
,{t,Range[0,1,0.01]}]
]


(* ::Subsection::Closed:: *)
(*hexToRGB*)


hexToRGB[hex_] := RGBColor @@ (IntegerDigits[ToExpression@StringReplace[hex, "#" -> "16^^"], 256, 3]/255.)


(* ::Subsection::Closed:: *)
(*monotonousDecreaseQ / monotonousIncreaseQ: Is a list monotonously decreasing / increasing?*)


monotonousIncreaseQ[list_]:=AllTrue[Differences[list],#>=0&]
monotonousDecreaseQ[list_]:=AllTrue[Differences[list],#<=0&]


(* ::Subsection::Closed:: *)
(*setRunTumbleTicks*)


(* ::Input:: *)
(*log10Q=True;*)
(*setRunTumbleTicks[log10Q];*)


Clear[setRunTumbleTicks];
setRunTumbleTicks[log10Q_:False]:=Module[{(*tickPositions,tickLabels,tickSize,*)ttTicks,trTicks},
	tickPositions=Flatten[Range[1.,9]#&/@{0.01,0.1,1,10}];
	tickSize=Table[If[Mod[Log10[val],1]==0,{0.016,0},{0.008,0}],{val,tickPositions}];
	tickLabels=Table[If[Mod[#,1]==0,Superscript[10,Round[#]],""]&@Log10[val],{val,tickPositions}];
	If[log10Q,tickPositions=Log10[tickPositions];];	
	
	ttTicks={tickPositions,tickLabels,tickSize}\[Transpose];
	trTicks={tickPositions,tickLabels,tickSize}\[Transpose];
	
	{trTicks,ttTicks}
]


(* ::Subsection::Closed:: *)
(*sortCleanTrajectories*)


(* ::Text:: *)
(*sort trajectories by decreasing length with minimum length of 10 s*)


(* ::Input:: *)
(*{trajectories,starts}=sortCleanTrajectories[trajectories0];*)


sortCleanTrajectories[trajectories0_]:=Module[{starts0,ends0,cleanTrajectories,minLengthSeconds,minLengthFrames,goodTrajsList,goodTrajectories,
goodStarts,sortByDecreasingLength,trajectories,starts},
	
	(** remove padding due to data storage scheme **)
	{starts0,ends0}=Transpose[
		{FirstPosition[#,x_/;x>0][[1]],
		Length[trajectories0]-FirstPosition[Reverse@#,x_/;x>0][[1]]+1
	}&/@(trajectories0\[Transpose][[All,All,1]])];
	cleanTrajectories=MapThread[#1[[#2;;#3]]&,{trajectories0\[Transpose],starts0,ends0}];
	
	(** get trajectories of minimum length **)
	minLengthSeconds=10.0;
	minLengthFrames=Round[minLengthSeconds/p["experiment","dt"]];
	goodTrajsList=(Length[#]>minLengthFrames)&/@cleanTrajectories;
	{goodTrajectories,goodStarts}=Pick[#,goodTrajsList]&/@{cleanTrajectories,starts0};
	
	(** sort trajectories by decreasing length **)
	sortByDecreasingLength=Reverse@Ordering[Length/@goodTrajectories];
	trajectories=goodTrajectories[[sortByDecreasingLength]];
	starts=goodStarts[[sortByDecreasingLength]];
	
	(** output **)
	{trajectories,starts}
]


(* ::Chapter:: *)
(*Figure 1*)


(* ::Section::Closed:: *)
(*Functions*)


(* ::Subsection::Closed:: *)
(*fig1*)


(* ::Input:: *)
(*fig1[p];*)


(* ::Subsubsection::Closed:: *)
(*Choose specific time trace interactively*)


(* ::Input:: *)
(*{dayIndexYoung,dayIndexOld}={1,6};*)
(*{dropletIndexYoung,dropletIndexOld}={3,9}(*{16,9}*);*)
(*{stepYoungStart,stepYoungEnd}=Round[{10,70}/p["experiment","dt"]];*)
(**)
(*dayIndices=Range[6];*)
(*allIntensities=getIntensitiesFromExperiment[p,dayIndices];*)
(*iTicks={0,1};*)
(*iTicks2={#,""}&/@iTicks;*)
(*tTicks=Table[t,{t,0,60,10}];*)
(*tTicks2={#,""}&/@tTicks;*)
(**)
(*Manipulate[*)
(*imax=Max[allIntensities[[dayIndexYoung,i]]];*)
(*trajectoryYoung={Range[0,Length[#]-1]p["experiment","dt"],#/imax}\[Transpose]&@allIntensities[[dayIndexYoung,i,stepYoungStart;;stepYoungEnd]];trajectoryOld={Range[0,Length[#]-1]p["experiment","dt"],#/imax}\[Transpose]&@allIntensities[[dayIndexOld,j]];ListLinePlot[{trajectoryYoung,trajectoryOld}*)
(*,Evaluate[p["plot","ops"]]*)
(*,AspectRatio->1/3*)
(*,Frame->True*)
(*,FrameLabel->{"time \!\(\*StyleBox[\"t\",\nFontSlant->\"Italic\"]\) (s)","Normalized intensity \!\(\*StyleBox[\"I\",\nFontSlant->\"Italic\"]\)"}*)
(*,FrameTicks->{{iTicks,iTicks2},{tTicks,tTicks2}}*)
(*,PlotRange->{{0,60},{0,1(*0.5*)}}*)
(*,PlotStyle->p["color","youngOld"]*)
(*]*)
(*,{{i,dropletIndexYoung},1,Length[allIntensities[[dayIndexYoung]]],1,Appearance->"Open"}*)
(*,{{j,dropletIndexOld},1,Length[allIntensities[[dayIndexOld]]],1,Appearance->"Open"}*)
(*,TrackedSymbols:>{i,j}*)
(*]*)


(* ::Subsubsection::Closed:: *)
(*Fig. 1d*)


(* ::Input:: *)
(*{trajectoryYoung,nMarkers,markerColors,markerSteps}=fig1PrepareData[p];*)
(*plot0=fig1dPlot[trajectoryYoung,markerSteps,markerColors,p];*)
(*plot=Show[plot0,ImageSize->90p["plot","mm"],ImageMargins->0(*,ImagePadding->{{11.5,2},{10,1}}mm*)]*)


(* ::Input:: *)
(*pthOutFig1=FileNameJoin[{p["path","figs"],"figure 1"}];*)
(*If[!DirectoryQ@#,CreateDirectory@#]&@pthOutFig1;*)
(*Export[FileNameJoin[{pthOutFig1,"fig1d.png"}],plot,ImageResolution->300]*)
(*Export[FileNameJoin[{pthOutFig1,"fig1d.pdf"}],plot]*)


(* ::Subsubsection::Closed:: *)
(*Fig. 1e*)


(* ::Input:: *)
(*{dayIndex,dropletIndex}={1,2};*)
(*{timeStart,timeDuration}={10,70};*)
(*singleDropletImages=getDropletFramesFromMovie[nMarkers,p,{timeStart,timeDuration,dayIndex,dropletIndex}]*)


(* ::Input:: *)
(*Export[FileNameJoin[{pthOutFig1,"fig1e.png"}],Grid[{singleDropletImages},Spacings->{0.1,0}]]*)


(* ::Subsubsection::Closed:: *)
(*code*)


fig1[p_]:=Module[{pthOutFig1,timeStart,timeEnd,trajectoryYoung,nMarkers,markerColors,markerSteps,plot0,plot1d,dayIndex,dropletIndex,timeStartMarkers,timeDurationMarkers,singleDropletImages},
	
	(** output directory **)
	pthOutFig1=FileNameJoin[{p["path","figs"],"figure 1"}];
	If[!DirectoryQ@#,CreateDirectory@#]&@pthOutFig1;
	
	{dayIndex,dropletIndex}={1,(*3*)2};   (** See later subsection: Choose specific time trace interactively **)
	{timeStart,timeEnd}={10,70};          (** in s **)
	{trajectoryYoung,nMarkers,markerColors,markerSteps}=fig1PrepareData[timeStart,timeEnd,dayIndex,dropletIndex,p];
	
	(** Figure 1d **)
	plot0=fig1dPlot[trajectoryYoung,markerSteps,markerColors,p];
	plot1d=Show[plot0,ImageSize->90p["plot","mm"],ImageMargins->0(*,ImagePadding->{{11.5,2},{10,1}}mm*)];
	Print[plot1d];
	Export[FileNameJoin[{pthOutFig1,"fig1d.png"}],plot1d,ImageResolution->300];
	Export[FileNameJoin[{pthOutFig1,"fig1d.pdf"}],plot1d];
	
	(** Figure 1e **)
	{timeStartMarkers,timeDurationMarkers}=timeStart+{0.1,0.9}(timeEnd-timeStart); 
	singleDropletImages=getDropletFramesFromMovie[nMarkers,p,{timeStartMarkers,timeDurationMarkers,dayIndex,dropletIndex}];
	Print[singleDropletImages];
	Export[FileNameJoin[{pthOutFig1,"fig1e.png"}],Grid[{singleDropletImages},Spacings->{0.1,0}]];
]


(* ::Subsection::Closed:: *)
(*fig1PrepareData*)


(* ::Input:: *)
(*{trajectoryYoung,nMarkers,markerColors,markerSteps}=fig1PrepareData[dayIndex,dropletIndex,p];*)


fig1PrepareData[timeStart_,timeEnd_,dayIndex_,dropletIndex_,p_]:=Module[{dropletIntensities,stepStart,stepEnd,imax,intensities,nMarkers,markerColors,markerSteps,timeIntensity},
	
	dropletIntensities=getIntensitySingleFromExperiment[p,dayIndex][[dropletIndex]];
	
	{stepStart,stepEnd}=Round[{timeStart,timeEnd}/p["experiment","dt"]];
	intensities=dropletIntensities[[stepStart;;stepEnd]];
	imax=1.2 Max[intensities]; (** factor 1.2 looks good for droplet 1:3 **)
	
	(** normalized intensity time trace **)
	timeIntensity={Range[0,Length[#]-1]p["experiment","dt"], #/imax}\[Transpose]&@intensities;
	
	(** markers **)
	nMarkers=4;
	markerColors=timeColormap/@Subdivide[0,1,nMarkers-1];
	markerSteps=Round[Subdivide[0.1,0.9,nMarkers-1]Length[timeIntensity]];   (** in frames after trajectory start **)
	
	(** output **)
	{timeIntensity,nMarkers,markerColors,markerSteps}
]


(* ::Subsection::Closed:: *)
(*fig1dPlot*)


(* ::Input:: *)
(*plot=fig1dPlot[trajectoryYoung,markerSteps,markerColors,p]*)


(* ::Subsubsection:: *)
(*code*)


Clear[fig1dPlot];
fig1dPlot[trajectoryYoung_,markerSteps_,markerColors_,p_]:=Module[{iTicks,iTicks2,tTicks,tTicks2},
	iTicks={0,1};
	iTicks2={#,""}&/@iTicks;
	tTicks=Table[t,{t,0,60,10}];
	tTicks2={#,""}&/@tTicks;
	
	(** output: plot **)
	ListLinePlot[{trajectoryYoung}
		,Evaluate[p["plot","ops"]]
		,AspectRatio->1/3
		,FrameTicks->{{iTicks,iTicks2},{tTicks,tTicks2}}
		,ImagePadding->{{34,12},{35,8}}
		,PlotRange->{{0,60},{0,1}}
		,PlotRangeClipping->False
		,PlotStyle->p["color","young"]
		,Epilog->{
			(** frame labels **)
			Text[Style[p["plot","timeLabel"],p["plot","fontSize"]],Scaled[{0.5,-0.27}]]
			,Rotate[Text[Style[p["plot","intensityLabel"],p["plot","fontSize"]],Scaled[{-0.08,0.5}]],90\[Degree]]
			(** color markers **)
			,{AbsolutePointSize[8],MapThread[{#1,Point[trajectoryYoung[[#2]]]}&,{markerColors,markerSteps}]}
		}
	]
]


(* ::Section:: *)
(*Plotting*)


(* ::Input:: *)
(*fig1[p];*)


(* ::Chapter:: *)
(*Figure 2*)


(* ::Section:: *)
(*Functions*)


(* ::Subsection:: *)
(*fig2*)


(* ::Input:: *)
(*fig2[p];*)


(* ::Subsubsection::Closed:: *)
(*debug: prepare data*)


(* ::Input:: *)
(*{timeStart,timeDuration}={10,30};*)
(*dayIndices={1,6};*)
(*dropletIndices={2,1,1,1,1,3}; (** found via Manipulate[] tool for Fig. 1d **)*)
(**)
(*{allIntensities,freqs,meanSpectra,stdSpectra,iAmplitudes,iWidths,allSpectra,pthOutputFig2}=fig2PrepareData[p];*)


(* ::Subsubsection::Closed:: *)
(*debug: Fig. 2a: Single timeseries for days {1, 6}*)


(* ::Input:: *)
(*timeseriesOverviewPlot=fig2aPlot[dayIndices,allIntensities,p,{timeStart,timeDuration,dropletIndices},markerColors,markerSteps];*)
(*plotYoungOld=Show[timeseriesOverviewPlot,ImageSize->90p["plot","mm"],ImagePadding->{{11.5,2},{11.5,1}}p["plot","mm"],ImageMargins->0]*)


(* ::Input:: *)
(*Export[FileNameJoin[{pthOutputFig2,"plot_timeseriesYoungOld.pdf"}],plotYoungOld]*)
(*Export[FileNameJoin[{pthOutputFig2,"plot_timeseriesYoungOld.png"}],plotYoungOld,ImageResolution->300]*)


(* ::Subsubsection::Closed:: *)
(*debug: Fig. 2a insets - single droplet images*)


(* ::Input:: *)
(*nSingleDropletImages=12;*)
(*{timeStartMarkers,timeDurationMarkers}=timeStart+{0.1,0.9}timeDuration; *)
(*singleDropletImagesAll=Table[*)
(*Print["dayIndex: ",dayIndex];singleDropletImages=getDropletFramesFromMovie[nSingleDropletImages,p,{timeStartMarkers,timeDurationMarkers,dayIndex,dropletIndices[[dayIndex]]}];Print[Grid[{#},Spacings->{0.1,0}]&@singleDropletImages];*)
(**)
(*(** output **)*)
(*singleDropletImages*)
(*,{dayIndex,dayIndices}];*)


(* ::Input:: *)
(*MapThread[Export[FileNameJoin[{pthOutputFig2,#1<>".png"}],Grid[{#2},Spacings->{0.1,0}]]&,{"dropletSnapshots_"<>ToString[#]<>"_"<>ToString[nSingleDropletImages]&/@dayIndices,singleDropletImagesAll}]//TableForm*)


(* ::Subsubsection::Closed:: *)
(*debug: Fig. 2b,c: Intensity spectrum for days {1,6}*)


(* ::Input:: *)
(*{plotStdQ,plotIndiQ}={False,False};*)
(*{plots0,stdplots0,indiPlots0}=fig2bcPlot[freqs,meanSpectra,iAmplitudes,iWidths,allSpectra,stdSpectra,p,{plotStdQ,plotIndiQ}];*)
(*plots=Show[#,ImageSize->44.5p["plot","mm"],ImageMargins->0]&/@plots0;*)
(*Grid[{plots}]*)
(**)
(*If[plotStdQ,*)
(*	Print["Plot with standard deviation:"];*)
(*	stdplots=Show[#,ImageSize->44.5p["plot","mm"],ImageMargins->0]&/@stdplots0;*)
(*	Grid[{stdplots}]*)
(*]*)
(**)
(*If[plotIndiQ,*)
(*	Print["Plot with all individual spectra:"];*)
(*	indiPlots=Show[#,ImageSize->44.5p["plot","mm"],ImageMargins->0]&/@indiPlots0;*)
(*	Grid[{indiPlots}]*)
(*]*)


(* ::Input:: *)
(*MapThread[*)
(*	Export[FileNameJoin[{pthOutputFig2,#2<>".png"}],#1,ImageResolution->300]*)
(*&,{plots,{"plot2_young_spectrum","plot2_old_spectrum"}}];*)
(*MapThread[*)
(*	Export[FileNameJoin[{pthOutputFig2,#2<>".pdf"}],#1]*)
(*&,{plots,{"plot2_young_spectrum","plot2_old_spectrum"}}];*)


(* ::Input:: *)
(*MapThread[*)
(*	Export[FileNameJoin[{pthOutputFig2,#2<>".png"}],#1,ImageResolution->300]*)
(*&,{stdplots,{"plot2_young_spectrum_std","plot2_old_spectrum_std"}}];*)


(* ::Input:: *)
(*MapThread[*)
(*	Export[FileNameJoin[{pthOutputFig2,#2<>".png"}],#1,ImageResolution->300]*)
(*&,{indiPlots,{"plot2_young_spectrum_indi","plot2_old_spectrum_indi"}}];*)


(* ::Subsubsection::Closed:: *)
(*code*)


fig2[p_]:=Module[{timeStart,timeDuration,dayIndices,dropletIndices,allIntensities,freqs,meanSpectra,stdSpectra,iAmplitudes,iWidths,allSpectra,pthOutputFig2,markerColors,markerSteps,
timeseriesOverviewPlot,fig2a,nSingleDropletImages,timeStartMarkers,timeDurationMarkers,singleDropletImagesAll,singleDropletImages,plotStdQ,plotIndiQ,plots0,stdplots0,indiPlots0,
figs2bc,stdplots,indiPlots},
	
	{timeStart,timeDuration}={10,30};
	dayIndices={1,6};
	dropletIndices={2,1,1,1,1,3}; (** found via Manipulate[] tool for Fig. 1d **)

	{allIntensities,freqs,meanSpectra,stdSpectra,iAmplitudes,iWidths,allSpectra,pthOutputFig2,markerColors,markerSteps}=fig2PrepareData[timeDuration,p];
	
	(** Fig. 2a **)
	timeseriesOverviewPlot=fig2aPlot[dayIndices,allIntensities,p,{timeStart,timeDuration,dropletIndices},markerColors,markerSteps];
	fig2a=Show[timeseriesOverviewPlot,ImageSize->90p["plot","mm"],ImagePadding->{{11.5,2},{11.5,1}}p["plot","mm"],ImageMargins->0];
	Print[fig2a];
	Export[FileNameJoin[{pthOutputFig2,"plot_timeseriesYoungOld.pdf"}],fig2a];
	Export[FileNameJoin[{pthOutputFig2,"plot_timeseriesYoungOld.png"}],fig2a,ImageResolution->300];
	
	(** Fig. 2a insets **)
	nSingleDropletImages=13;
	{timeStartMarkers,timeDurationMarkers}=timeStart+{0.1,0.9}timeDuration; 
	singleDropletImagesAll=Table[
		Print["dayIndex: ",dayIndex];
		singleDropletImages=getDropletFramesFromMovie[nSingleDropletImages,p,{timeStartMarkers,timeDurationMarkers,dayIndex,dropletIndices[[dayIndex]]}];
		Print[Grid[{#},Spacings->{0.1,0}]&@singleDropletImages];
	
		(** output **)
		singleDropletImages
	,{dayIndex,dayIndices}];
	MapThread[Export[FileNameJoin[{pthOutputFig2,#1<>".png"}],Grid[{#2},Spacings->{0.1,0}]]&,{"dropletSnapshots_"<>ToString[#]<>"_"<>ToString[nSingleDropletImages]&/@dayIndices,singleDropletImagesAll}];
	
	(** Figs. 2b,c **)
	{plotStdQ,plotIndiQ}={False,False};
	{plots0,stdplots0,indiPlots0}=fig2bcPlot[freqs,meanSpectra,iAmplitudes,iWidths,allSpectra,stdSpectra,p,{plotStdQ,plotIndiQ}];
	figs2bc=Show[#,ImageSize->44.5p["plot","mm"],ImageMargins->0]&/@plots0;
	Print[Grid[{figs2bc}]];
	MapThread[
		Export[FileNameJoin[{pthOutputFig2,"plot2_"<>#2<>"_spectrum.png"}],#1,ImageResolution->300]
	&,{figs2bc,{"young","old"}}];
	MapThread[
		Export[FileNameJoin[{pthOutputFig2,"plot2_"<>#2<>"_spectrum.pdf"}],#1]
	&,{figs2bc,{"young","old"}}];
	
	If[plotStdQ,
		Print["Plot with standard deviation:"];
		stdplots=Show[#,ImageSize->44.5p["plot","mm"],ImageMargins->0]&/@stdplots0;
		Print[Grid[{stdplots}]];
		MapThread[
			Export[FileNameJoin[{pthOutputFig2,#2<>".png"}],#1,ImageResolution->300]
		&,{stdplots,{"plot2_young_spectrum_std","plot2_old_spectrum_std"}}];
	];
	
	If[plotIndiQ,
		Print["Plot with all individual spectra:"];
		indiPlots=Show[#,ImageSize->44.5p["plot","mm"],ImageMargins->0]&/@indiPlots0;
		Print[Grid[{indiPlots}]];
		MapThread[
			Export[FileNameJoin[{pthOutputFig2,#2<>".png"}],#1,ImageResolution->300]
		&,{indiPlots,{"plot2_young_spectrum_indi","plot2_old_spectrum_indi"}}];
	];
]


(* ::Subsection::Closed:: *)
(*fig2PrepareData: get intensity time series and average spectrum for each day*)


(* ::Input:: *)
(*{allIntensities,freqs,meanSpectra,stdSpectra,iAmplitudes,iWidths,allSpectra,pthOutputFig2,markerColors,markerSteps}=fig2PrepareData[timeDuration,p]*)


(* ::Subsubsection:: *)
(*code*)


Clear[fig2PrepareData];
fig2PrepareData[timeDuration_,p_]:=Module[
{allIntensities,freqs,weightedMeanSpectra,weightedStdSpectra,allWeightedSpectra,fMaxIndex,indexStart,indexEnd,iFits,iAmplitudes,iWidths,powerSpectralDensities,pthOutputFig2,
nMarkers,markerColors,markerSteps},

	{allIntensities,freqs,weightedMeanSpectra,weightedStdSpectra,allWeightedSpectra}=calculateDailyMeanSpectra[p][[{1,2,6,7,8}]];
	
	(** fit mean power spectral densities **)
	fMaxIndex=FirstPosition[freqs,_?(#>=p["fourier","fMax"]&)][[1]];
	{indexStart,indexEnd}={1,fMaxIndex};
	(** Lorentz with peak at 0 Hz **)
	iFits=NonlinearModelFit[{freqs,#/Max[#]}\[Transpose][[indexStart;;indexEnd]],{aa/(1+(\[Omega]/bb)^2),aa>0,bb>0},{{aa,1.0},{bb,0.1}},\[Omega],MaxIterations->1000]&/@weightedMeanSpectra;   
	iAmplitudes=(#["BestFitParameters"][[All,2]]&/@iFits)\[Transpose][[1]];
	iWidths=(#["BestFitParameters"][[All,2]]&/@iFits)\[Transpose][[2]];
	
	(** output directory **)
	pthOutputFig2=FileNameJoin[{p["path","figs"],"figure 2"}];
	If[!DirectoryQ@#,CreateDirectory@#]&@pthOutputFig2;
	
	(** markers **)
	nMarkers=4;
	markerColors=timeColormap/@Subdivide[0,1,nMarkers-1];
	markerSteps=Round[Subdivide[0.1,0.9,nMarkers-1]timeDuration/p["experiment","dt"]];   (** in frames after trajectory start **)
	
	(** output **)
	{allIntensities,freqs,weightedMeanSpectra,weightedStdSpectra,iAmplitudes,iWidths,allWeightedSpectra,pthOutputFig2,markerColors,markerSteps}
]


(* ::Subsection::Closed:: *)
(*fig2aPlot*)


(* ::Input:: *)
(*dayIndices={1,6};*)
(*timeseriesOverviewPlots=fig2aPlot[dayIndices,allIntensities,p,{timeStart,timeDuration,dropletIndices},markerColors,markerSteps];*)


(* ::Subsubsection::Closed:: *)
(*debug: font family for serif capital I*)


(* ::Input:: *)
(*Manipulate[*)
(*Row[{"Intensity ",Style["I",FontFamily->font,FontSlant->Italic]}]*)
(*,{font,$FontFamilies}]*)


(* ::Subsubsection:: *)
(*code*)


Clear[fig2aPlot];
fig2aPlot[dayIndices_,allIntensities_,p_,{timeStart_,timeDuration_,dropletIndices_},markerColors_,markerSteps_]:=Module[{stepStart,stepEnd,imin,imax,
tTicks,tTicks2,dropletIndex,intensity,time,space,xLabel,yLabel,timeIntensities},
	
	(** options **)
	{stepStart,stepEnd}=Round[{timeStart,timeStart+timeDuration}/p["experiment","dt"]];
	{imin,imax}={0.0,0.35};   (** y plotrange **)
	tTicks=Range[0,30,10];
	tTicks2={#,""}&/@tTicks;
	
	(** plot labels **)
	space=-8;
	xLabel=Framed[p["plot","timeLabel"],FrameStyle->None,ContentPadding->False,FrameMargins->{{0,0},{0,space}}];
	yLabel=p["plot","intensityLabel"];
	
	(** timeseries plot **)
	timeIntensities=Table[
		If[!MemberQ[{1,6},dayIndex],Print["Warning: This code is only tested for dayIndices \[Element] {1,6}"];];
		dropletIndex=dropletIndices[[dayIndex]];
	
		(** create time series plots **)
		intensity=Rescale[allIntensities[[dayIndex,dropletIndex,stepStart;;stepEnd]],{imin,imax}];
		time=Range[0,Length[#]-1]p["experiment","dt"]&@intensity;
		
		(** output **)
		{time,intensity}\[Transpose]
	,{dayIndex,dayIndices}];
	
	(** output **)
	ListLinePlot[timeIntensities
		,Evaluate[p["plot","ops"]]
		,AspectRatio->1/2.5
		,FrameLabel->{xLabel,yLabel}
		,FrameTicks->{{None,None},{tTicks,tTicks2}}
		(*,ImagePadding\[Rule]{{35,11},{15,2}}*)
		,PlotRange->{{0,timeDuration},{0,1}}
		,PlotRangePadding->{{Scaled[0.05],Scaled[0.05]},{Scaled[0.15],Scaled[0.15]}}
		,PlotStyle->Evaluate[Directive[p["color","ages"][[#]],AbsoluteThickness[p["plot","lineThickness"]]]&/@dayIndices]
		,Epilog->{
			(** frame labels **)
(*				Text[Style[p["plot","timeLabel"],p["plot","fontSize"]],Scaled[{0.5,-0.27}]]
			,Rotate[Text[Style[p["plot","intensityLabel"],p["plot","fontSize"]],Scaled[{-0.08,0.5}]],90\[Degree]]*)
			(** color markers **)
			(*,*){AbsolutePointSize[4],Table[MapThread[{#1,Point[timeIntensity[[#2]]]}&,{markerColors,markerSteps}],{timeIntensity,timeIntensities}]}
		}
	]
]


(* ::Subsection::Closed:: *)
(*fig2bcPlot*)


(* ::Input:: *)
(*{plotStdQ,plotIndiQ}={True,False};*)
(*{plots,stdplots,indiPlots}=fig2bcPlot[freqs,meanSpectra,iAmplitudes,iWidths,allSpectra,stdSpectra,p,{plotStdQ,plotIndiQ}];*)


(* ::Subsubsection::Closed:: *)
(*view individual trajectories and power spectra*)


(* ::Input:: *)
(*dayIndex=1;*)
(*Manipulate[*)
(*	label=Style["day: "<>ToString[IntegerString[dayIndex]]<>", droplet: "<>ToString[IntegerString[dropletIndex]],Black];*)
(*	p1=ListLinePlot[allIntensities[[dayIndex,dropletIndex]],(*PlotRange->{{-0.01 fMax,1.01fMax},{-0.01,1.05}},*)Frame->True,AspectRatio->1,FrameLabel->{"time \!\(\**)
(*	StyleBox[\"t\",\nFontSlant->\"Italic\"]\) (s)","Intensity \!\(\**)
(*	StyleBox[\"I\",\nFontSlant->\"Italic\"]\)"},FrameStyle->Directive[Black,18,AbsoluteThickness[2]],ImageSize->{Automatic,300}];*)
(*	p2=ListLinePlot[{freqs,#/0.05}\[Transpose](*\[LeftDoubleBracket]2;;\[RightDoubleBracket]*)&@allSpectra[[dayIndex,dropletIndex]],PlotRange->{{-0.01 p["fourier","fMax"],1.01p["fourier","fMax"]},{-0.01,1.05}},Frame->True,AspectRatio->1,FrameLabel->{"frequency \!\(\**)
(*	StyleBox[\"f\",\nFontSlant->\"Italic\"]\) (Hz)","norm. PSD (a.u.)"},PlotStyle->(ColorData[{"BlueGreenYellow","Reverse"}][#]&/@Subdivide[Length[meanSpectra]-1]),FrameStyle->Directive[Black,18,AbsoluteThickness[2]],ImageSize->{Automatic,300}];*)
(*	Grid[{{label,SpanFromLeft},{p1,p2}}]*)
(*,{dropletIndex,1,Length[allSpectra[[dayIndex]]],1,Appearance->"Open"}]*)


(* ::Subsubsection::Closed:: *)
(*daily mean spectra*)


(* ::Input:: *)
(*filterQ=True;*)
(*plot=ListLinePlot[{freqs,#/Max[#]}\[Transpose](*\[LeftDoubleBracket]2;;\[RightDoubleBracket]*)&/@meanSpectra,PlotRange->{{0,p["fourier","fMax"]},{0,1}},PlotLegends->Table["day "<>ToString[i],{i,6}],Frame->True,AspectRatio->1,FrameLabel->{"Frequency \!\(\**)
(*StyleBox[\"f\",\nFontSlant->\"Italic\"]\) (Hz)","Intensity\nspectrum"},PlotStyle->(ColorData[{"BlueGreenYellow","Reverse"}][#]&/@Subdivide[Length[meanSpectra]-1]),FrameStyle->Directive[Black,18,AbsoluteThickness[2]],PlotLabel->Style["Individually normalized daily mean spectra",Black]]*)
(*Export[FileNameJoin[{pthOutput,"meanSpectraAll"<>If[filterQ,"_filtered",""]<>"."<>#}],plot]&/@{"png","pdf"}//TableForm*)


(* ::Input:: *)
(*Table[*)
(*ListLinePlot[{freqs,#/Max[#]}\[Transpose](*\[LeftDoubleBracket]2;;\[RightDoubleBracket]*)&@meanSpectra[[dayIndex]]*)
(*,AspectRatio->1*)
(*,Frame->True*)
(*,FrameLabel->{"Frequency \!\(\**)
(*StyleBox[\"f\",\nFontSlant->\"Italic\"]\) (Hz)","Intensity\nspectrum"}*)
(*,FrameStyle->Directive[Black,18,AbsoluteThickness[2]]*)
(*,ImageSize->{Automatic,300}*)
(*,PlotRange->{{-0.01,1.01}p["fourier","fMax"],{-0.02,1.05}}*)
(*,PlotStyle->(ColorData[{"BlueGreenYellow","Reverse"}][#]&@((dayIndex-1)/(p["experiment","nDays"]-1)))*)
(*]*)
(*,{dayIndex,Range[p["experiment","nDays"]]}]*)


(* ::Subsubsection::Closed:: *)
(*code*)


fig2bcPlot[freqs_,meanSpectra_,iAmplitudes_,iWidths_,allSpectra_,stdSpectra_,p_,{plotStdQ_,plotIndiQ_}]:=Module[{fTicks,fTicks2,pTicks,pTicks2,yMin,yMax,space,xLabel,yLabel,ops,dayIndices,dayIndex,plots
,stdplots,indiPlots,color,meanSpectrum,fMaxIndex,iAmplitude,iWidth,fwhmPoint,psdMax,meanPSD,pMean,meanPSDpSTD,meanPSDmSTD,pMeanStd,pMeanAndIndividuals,imageSize},

	(** tick options **)
	fTicks=Table[If[Mod[t,0.5]==0,{t,DecimalForm[t,{2,1}],{0.03,0}},{t,"",{0.015,0}}],{t,0,1,0.1}];
	fTicks2={#[[1]],"",#[[3]]}&/@fTicks;
	pTicks=Table[{t,DecimalForm[t,{2,1}]},{t,0,1,0.5}];
	pTicks2={#[[1]],""}&/@pTicks;
	
	(** plot range **)
	{yMin,yMax}={-0.02,4.5};
	
	(** plot labels **)
	space=-8;
	xLabel=Framed[p["plot","frequencyLabel"],FrameStyle->None,ContentPadding->False,FrameMargins->{{0,0},{0,space}}];
	yLabel=p["plot","spectrumLabel"];
	
	
	(** plot options **)
	ops={
		Evaluate[p["plot","ops"]]
		,AspectRatio->1
		,FrameLabel->{xLabel,yLabel}
		,FrameTicks->{{pTicks,pTicks2},{fTicks,fTicks2}}
		(*,PlotRange->{{-0.01,1.01}fMax,{yMin,yMax}}*)
	};
	
	
	(** select days to plot **)
	dayIndices={1,6};
	
	(** plots for both days **)
	{plots,stdplots,indiPlots}=Table[
		color=p["color","youngOld"][[plotIndex]];
		dayIndex=dayIndices[[plotIndex]];
		meanSpectrum=meanSpectra[[dayIndex]];

		(** extract data from fitted spectra **)
		{iAmplitude,iWidth}={iAmplitudes[[dayIndex]],iWidths[[dayIndex]]};
		fwhmPoint={iWidth,0.5};

		psdMax=iAmplitude Max[meanSpectrum];
		meanPSD={freqs,#/psdMax}\[Transpose]&@meanSpectrum;
		pMean=Show[
			ListLinePlot[meanPSD
				,PlotRange->{{-0.01,1.01}p["fourier","fMax"],All}
				,Evaluate[ops]
				,PlotStyle->color
				,Epilog->{
					Red,AbsolutePointSize[7.5],Point[#],AbsoluteThickness[p["plot","lineThickness"]],Dashed
					,HalfLine[{#,{0,0.5}}],HalfLine[{#,{#[[1]],0}}]
					,Text[Style["\!\(\*StyleBox[SubscriptBox[\"f\", \"HM\"],\nFontSlant->\"Italic\"]\)",p["plot","fontSize"],p["plot","fontFamily"],Red],#,{-1.6,-1}]}&@fwhmPoint
			]
		];
		
		If[plotStdQ,
			meanPSDpSTD={freqs,(#[[1]]+#[[2]])/Max[#[[1]]]}\[Transpose]&@{meanSpectrum,stdSpectra[[dayIndex]]};
			meanPSDmSTD={freqs,(#[[1]]-#[[2]])/Max[#[[1]]]}\[Transpose]&@{meanSpectrum,stdSpectra[[dayIndex]]};
			meanPSDmSTD[[All,2]]=Max[#,0]&/@meanPSDmSTD[[All,2]];
			pMeanStd=Show[
				ListLinePlot[{meanPSD,meanPSDpSTD,meanPSDmSTD}
					,PlotRange->{{-0.01,1.01}p["fourier","fMax"],{yMin,yMax}}
					,Evaluate[ops]
					,PlotStyle->{color,None,None}
					,Filling->{2->{3}}
					,FillingStyle->Lighter@color
					,Epilog->{Red,AbsolutePointSize[7.5],Point[#],AbsoluteThickness[p["plot","lineThickness"]],Dashed,HalfLine[{#,{0,0.5}}],HalfLine[{#,{#[[1]],0}}]}&@fwhmPoint
				]
			]
		];
		
		If[plotIndiQ,
			pMeanAndIndividuals=Show[
				ListLinePlot[{freqs,#/psdMax}\[Transpose],PlotRange->{{-0.01,1.01}p["fourier","fMax"],{yMin,yMax}},Evaluate[ops],PlotStyle->LightGray]&/@allSpectra[[dayIndex]]
				,ListLinePlot[meanPSD,PlotStyle->color,Evaluate[ops]]
				,Epilog->{Red,AbsolutePointSize[7.5],Point[#],AbsoluteThickness[p["plot","lineThickness"]],Dashed,HalfLine[{#,{0,0.5}}],HalfLine[{#,{#[[1]],0}}]}&@fwhmPoint
			]
		];
		
		(** output **)
		{pMean,pMeanStd,pMeanAndIndividuals}
	,{plotIndex,{1,2}}]\[Transpose];

	(** output **)
	{plots,stdplots,indiPlots}
]


(* ::Section:: *)
(*Plotting*)


(* ::Input:: *)
(*fig2[p];*)


(* ::Chapter:: *)
(*Figure 3*)


(* ::Section:: *)
(*Functions*)


(* ::Subsection::Closed:: *)
(*createMovieFramesSimulation*)


(* ::Subsubsection::Closed:: *)
(*preliminary*)


(* ::Text:: *)
(*requires droplet_optics code to be run in advance to compute liquid crystal alignments*)


(* ::Input:: *)
(*solution=solveHemisphereLC[];*)


(* ::Subsubsection::Closed:: *)
(*test: few droplet frames*)


(* ::Input:: *)
(*{tmax,tfinalSim}={0.1,200};*)
(*dt=0.01;*)
(*{timeDuration,tStart,stepStart,stepEnd}=fig3SetTime[dt,tmax];*)
(**)
(*trajectory=trajectories[[trajIndex,stepsTransient;;stepsTransient+Round[tmax/tfinalSim nSteps];;1]];*)
(*ageString="young";*)
(**)
(*ageColor=p["color",ageString];*)
(*trajDataAll=trajDataAllYoung;*)


(* ::Input:: *)
(*{runColor,tumbleColor}={#,Darker[#,0.5]}&@ageColor;*)
(**)
(*counter=0;*)
(*outputDirectory=FileNameJoin[{NotebookDirectory[],"mmamovie_runTumbleSimulation_"<>ageString}];*)
(*If[!DirectoryQ[#],CreateDirectory[#]]&@outputDirectory;*)
(**)
(*positions=trajectory[[All,1;;3]];*)
(*thetas=ArcCos/@positions[[All,3]];   (** 0 - top pole (z=+1), \[Pi] - bottom pole (z=-1) **)*)
(*phis=phaseUnwrap[ArcTan[#[[1]],#[[2]]]&/@positions[[All,{1,2}]]];*)
(*orientations=trajectory[[All,4;;6]];*)
(*rtDatas=trajectory[[All,7]];*)
(**)
(*{timeDuration,tStart,stepStart,stepEnd}=fig3SetTime[dt,tmax];*)
(*nIntensityImages=stepEnd-stepStart+1;*)


(* ::Input:: *)
(*trajIndex=1;*)
(*{showTumbleShadingQ,colorizeTumblesQ}={True,True};*)
(*\[Theta]PlotMaxDeg=60;*)
(*\[Phi]OffsetTimeseries=0.0;*)
(**)
(*{timeDuration,tStart,stepStart,stepEnd}=fig3SetTime[dt];*)
(*stepPlot=100;*)
(*showTumbleShadingQ=False;*)
(**)
(**)
(*(** frame creation loop **)*)
(*Do[*)
(*	position=positions[[t]];*)
(*	orientation=orientations[[t]];*)
(*	theta=thetas[[t]];*)
(*	phi=phis[[t]];*)
(*	rtData=If[#<0.5,0,1]&/@rtDatas[[;;t]];*)
(*	{front,end}=Map[(1+rEcoliShort)position+# ((*rEcoliLong*)2rEcoliShort-rEcoliShort) orientation&,{1,-1}];*)
(*	*)
(*	(** liquid crystal droplet schematic subplot **)*)
(*	liquidCrystals=renderLCs[solution,{theta,phi}];*)
(**)
(*	droplet=ParametricPlot3D[{*)
(*			Evaluate[RotationMatrix[{{0,0,1},position}] . er[\[Theta],\[Phi]]]*)
(*			,Evaluate[RotationMatrix[{{0,0,1},position}] . er[\[Pi]-\[Theta],\[Phi]]]*)
(*		},{\[Theta],0,0.5\[Pi]},{\[Phi],0,2\[Pi]}*)
(*		,Mesh->None*)
(*		,Axes->False*)
(*		,Boxed->False*)
(*		,PlotPoints->50*)
(*		,PlotStyle->Evaluate[Directive[Opacity[0.2],Specularity[White,65],#]&/@{Lighter[Red],Lighter[Gray,0.7]}]*)
(*	];*)
(*	*)
(*	plotDropletLC=Show[droplet,liquidCrystals,PlotRange->All,Lighting->"Neutral"];*)
(*	(*Graphics3D[{Opacity[0.1],Sphere[]}]*)*)
(*	*)
(*	(** at origin **)*)
(*	(*plotAxesCross=Graphics3D[{Black,Arrow[Tube[{{0,0,0},1.3#}]]&/@{{1,0,0},{0,1,0},{0,0,1}}}];*)*)
(*	*)
(*	(** bottom left **)*)
(*	origin={0.8 1.0,0.8*-1.0,-1};*)
(*	plotAxesCross=Graphics3D[{Black,Arrowheads[0.03],Arrow[Tube[{origin,origin+0.35#}]]&/@{{-1,0,0},{0,1,0},{0,0,0.9 1}}}];*)
(**)
(*	plotBacteriumSphere=Show[*)
(*		plotDropletLC*)
(*		,plotAxesCross*)
(*		,Graphics3D[{ageColor,AbsoluteThickness[5],Line[#&/@trajectory[[;;t,1;;3]],VertexColors->(If[#==0,runColor,tumbleColor]&/@rtData)]}]	*)
(*		,Graphics3D[{If[rtData[[t]]==0,runColor,tumbleColor],CapsuleShape[{front,end},rEcoliShort]}]*)
(*		,Axes->False*)
(*		,Boxed->False*)
(*		,ImageSize->500*)
(*		,PlotRange->1.4{{-1,1},{-1,1},{-1,1}}*)
(*		,ViewVertical->{0.376,0.248,0.893}*)
(*		,ViewPoint->{2.813,1.728,0.742}*)
(*	];*)
(**)
(*	bacteriumAge=If[ageString=="young","Fresh E. coli","Starved E. coli"];*)
(*	timeLabel=Style[bacteriumAge<>"\ntime \!\(\*StyleBox[\"t\",\nFontSlant->\"Italic\"]\) = "<>ToString@DecimalForm[t dt,{5,Round[-Log10[dt]]}]<>"\[ThinSpace]s",FontFamily->"Helvetica",26];*)
(*	*)
(*	(** compose frame **)*)
(*	plot=plotBacteriumSphere;*)
(*	*)
(*	(** export frame **)*)
(*	Export[FileNameJoin[{outputDirectory,"p_"<>IntegerString[counter++,10,IntegerLength[Length[trajectory]]+1]<>".png"}],plot];*)
(**)
(*,{t,1,Length[trajectory],1}];*)


(* ::Subsubsection::Closed:: *)
(*test: few full frames*)


(* ::Input:: *)
(*{tmax,tfinalSim}={0.1,200};*)
(*dt=0.01;*)
(*{timeDuration,tStart,stepStart,stepEnd}=fig3SetTime[dt,tmax];*)
(**)
(*trajectory=trajectories[[trajIndex,stepsTransient;;stepsTransient+Round[tmax/tfinalSim nSteps];;1]];*)
(*ageString="young";*)
(**)
(*ageColor=p["color",ageString];*)
(*trajDataAll=trajDataAllYoung;*)


(* ::Input:: *)
(*{runColor,tumbleColor}={#,Darker[#,0.5]}&@ageColor;*)
(**)
(*counter=0;*)
(*outputDirectory=FileNameJoin[{NotebookDirectory[],"mmamovie_runTumbleSimulation_"<>ageString}];*)
(*If[!DirectoryQ[#],CreateDirectory[#]]&@outputDirectory;*)
(**)
(*positions=trajectory[[All,1;;3]];*)
(*thetas=ArcCos/@positions[[All,3]];   (** 0 - top pole (z=+1), \[Pi] - bottom pole (z=-1) **)*)
(*phis=phaseUnwrap[ArcTan[#[[1]],#[[2]]]&/@positions[[All,{1,2}]]];*)
(*orientations=trajectory[[All,4;;6]];*)
(*rtDatas=trajectory[[All,7]];*)
(**)
(*{timeDuration,tStart,stepStart,stepEnd}=fig3SetTime[dt,tmax];*)
(*nIntensityImages=stepEnd-stepStart+1;*)


(* ::Input:: *)
(*intensityImages0=getIntensityImages[thetas,phis,{stepStart,stepEnd},p,nIntensityImages];*)
(**)
(*trajIndex=1;*)
(*{showTumbleShadingQ,colorizeTumblesQ}={True,True};*)
(*\[Theta]PlotMaxDeg=60;*)
(*\[Phi]OffsetHemisphere=-0.5\[Pi];*)
(*\[Phi]OffsetTimeseries=0.0;*)
(**)
(*{timeDuration,tStart,stepStart,stepEnd}=fig3SetTime[dt];*)
(*stepPlot=100;*)
(*showTumbleShadingQ=False;*)
(**)
(*{showThetaMaxQ,sData,thetaData,phiData,iData,runColor,tumbleColor,ipx,ops,tTicks,tTicks2,\[Theta]maxModel,rtData,tumbleStarts,tumbleLengths}=fig3PrepareData[stepPlot,trajDataAll,trajIndex,ageString,colorizeTumblesQ,u0Dim,dt,tStart,stepStart,stepEnd,timeDuration,showTumbleShadingQ,p];*)


(* ::Subsubsection::Closed:: *)
(*debug: grid placement*)


(* ::Input:: *)
(*Print["movie frame generation took: ",AbsoluteTiming[*)
(*createMovieFramesSimulation[trajectories[[trajIndex,stepsTransient;;stepsTransient+Round[tmax/tfinalSim nSteps]]],"young",dt,trajDataAllYoung,timeDataYoung,freqDataYoung,pthFrames,p,tmax*)
(*];*)
(*][[1]],"\[ThinSpace]s"];*)
(*t=1;*)
(*Grid[{*)
(*	{timeLabel,SpanFromLeft},*)
(*    {mechanicsLabel,opticsLabel},*)
(*{*)
(*				Show[plotBacteriumSphere,ImageMargins->{{40,0},{0,0}}]*)
(*				,Rasterize[Show[intensityImages0[[t]],ImageSize->296,ImageMargins->{{80,0},{0,0}}]]*)
(*}*)
(*	,{*)
(*		Column[{pRT,pTheta,pPhi},Alignment->Center,Spacings->0,Frame->All]*)
(*		,Column[{pI,pSpectrum},Alignment->Center,Spacings->0,Frame->All]*)
(*	}*)
(*},Spacings->0,Frame->All]*)


(* ::Subsubsection::Closed:: *)
(*test: single frame*)


(* ::Input:: *)
(*createMovieFramesSimulation[*)
(*	trajectories[[trajIndex,stepsTransient;;stepsTransient+1]]*)
(*,"young",dt,trajDataAllYoung,timeDataYoung,freqDataYoung,pthFrames,tmax]*)


(* ::Subsubsection::Closed:: *)
(*code**)


(* ::Text:: *)
(*TODO: only a single trajectory data source, instead of two!*)


Clear[createMovieFramesSimulation];
createMovieFramesSimulation[trajectory_,ageString_,dt_,trajDataAll_,timeData_,freqData_,outputDirectory_,p_,tmax_:20,frameStep_:1]:=Module[{color,counter,rtData,position,theta,phi,
orientation,front,end,runColor,tumbleColor,bacteriumAge,timeLabel,plot,liquidCrystals,droplet,origin,plotAxesCross,plotDropletLC(*,plotBacteriumSphere*)
,(*pRT,pTheta,pPhi,*)opsMovie},
	
	counter=0;
	If[!DirectoryQ[#],CreateDirectory[#]]&@outputDirectory;
	
	(** trajectory source 1 **)
	positions=trajectory[[All,1;;3]];
	thetas=ArcCos/@positions[[All,3]];   (** 0 - top pole (z=+1), \[Pi] - bottom pole (z=-1) **)
	phis=phaseUnwrap[ArcTan[#[[1]],#[[2]]]&/@positions[[All,{1,2}]]];
	orientations=trajectory[[All,4;;6]];
	rtDatas=trajectory[[All,7]];
	
	{timeDuration,tStart,stepStart,stepEnd}=fig3SetTime[dt,tmax];
	(** TODO: account for frameStep option **)
	nIntensityImages=stepEnd-stepStart+1;
	Print["finished loading simulated intensity images: ",AbsoluteTiming[
		intensityImages0=getIntensityImages[thetas,phis,{stepStart,stepEnd},p,nIntensityImages];
	][[1]],"\[ThinSpace]s"];
	
	trajIndex=1;
	{showTumbleShadingQ,colorizeTumblesQ}={True,True};
	\[Theta]PlotMaxDeg=60;
	\[Phi]OffsetTimeseries=0.0;
	
	stepPlot=100;
	showTumbleShadingQ=False;
	
	(** trajectory source 2 **)
	{showThetaMaxQ,sData,thetaData,phiData,iData,runColor,tumbleColor,ipx,ops,tTicks,tTicks2,\[Theta]maxModel,rtData,tumbleStarts,tumbleLengths}=fig3PrepareData[stepPlot,trajDataAll,trajIndex,ageString,colorizeTumblesQ,u0Dim,dt,tStart,stepStart,stepEnd,timeDuration,showTumbleShadingQ,p];
	
	
	(** frame creation loop **)
	Do[
		position=positions[[t]];
		orientation=orientations[[t]];
		theta=thetas[[t]];
		phi=phis[[t]];
		rtData=If[#<0.5,0,1]&/@rtDatas[[;;t]];
		{front,end}=Map[(1+rEcoliShort)position+# ((*rEcoliLong*)2rEcoliShort-rEcoliShort) orientation&,{1,-1}];
		
		(** liquid crystal droplet schematic subplot **)
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
	
		plotBacteriumSpherePadded=Show[
			plotDropletLC
			,plotAxesCross
			,Graphics3D[{runColor,AbsoluteThickness[5],Line[#&/@trajectory[[;;t,1;;3]],VertexColors->(If[#==0,runColor,tumbleColor]&/@rtData)]}]	
			,Graphics3D[{If[rtData[[t]]==0,runColor,tumbleColor],CapsuleShape[{front,end},rEcoliShort]}]
			,Axes->False
			,Boxed->False
			,ImageSize->500
			,PlotRange->1.4{{-1,1},{-1,1},{-1,1}}
			,ViewVertical->{0.376,0.248,0.893}
			,ViewPoint->{2.813,1.728,0.742}
		];
		plotBacteriumSphere=ImageTake[plotBacteriumSpherePadded,{44,438},{165,553}];
	
		bacteriumAge=If[ageString=="young","Fresh","Starved"]<>" E. coli";
		timeLabel=Style[bacteriumAge<>"\ntime \!\(\*StyleBox[\"t\",\nFontSlant->\"Italic\"]\) = "<>ToString@DecimalForm[(t-1) dt,{5,Round[-Log10[dt]]}]<>"\[ThinSpace]s",FontFamily->"Helvetica",22,Black];
		mechanicsLabel=Style["    Bacterial kinematics",FontFamily->"Helvetica",18,Black];
		opticsLabel=Style["       Droplet optics",FontFamily->"Helvetica",18,Black];
		
		(** update plot option **)
		{ipx,ops}=fig3CommonPlotOptions[t,sData,dt,runColor,tumbleColor,p];
		
		(** plot: RT state **)
		showTumbleShadingQ=False;
		opsMovie={AspectRatio->1/3,PlotStyle->AbsolutePointSize[3],ImageSize->55p["plot","mm"](*,ImagePadding->{ipx,{5,8}}*)};
		pRT=fig3PlotRT[timeData,sData,stepStart,stepEnd,ipx,ops,tTicks2,showTumbleShadingQ,tumbleStarts,tumbleLengths,p,opsMovie];
		
		(** plots: polar and azimuthal angle **)
		opsMovie={
			AspectRatio->1/3
			,PlotStyle->AbsoluteThickness[3]
			,ImagePadding->{ipx,{11,14}}
			,ImageSize->55p["plot","mm"]
			,Epilog->{Rotate[Text[Style["\!\(\*StyleBox[\"\[Theta]\",\nFontSlant->\"Italic\"]\) (\[Degree])",p["plot","fontSize"]],Scaled[{-0.27,0.5}]],90\[Degree]]}
		};
		pTheta=fig3PlotTheta[\[Theta]PlotMaxDeg,timeData,thetaData,stepStart,stepEnd,ipx,ops,tTicks2,showTumbleShadingQ,tumbleStarts,tumbleLengths,showThetaMaxQ,\[Theta]maxModel,p,opsMovie];
		
		opsMovie={
			AspectRatio->1/3
			,ImagePadding->{ipx,{36,10}}
			,ImageSize->55p["plot","mm"]
			,Epilog->{
				Text[Style[p["plot","timeLabel"],p["plot","fontSize"]],Scaled[{0.5,-0.55}]]
				,Rotate[Text[Style["\!\(\*StyleBox[\"\[Phi]\",\nFontSlant->\"Italic\"]\) (\[Degree])",p["plot","fontSize"]],Scaled[{-0.27,0.5}]],90\[Degree]]
			}
			,PlotStyle->AbsoluteThickness[3]
		};
		pPhi=fig3PlotPhi[timeData,phiData,stepStart,stepEnd,\[Phi]OffsetTimeseries,ipx,ops,tTicks,tTicks2,showTumbleShadingQ,tumbleStarts,tumbleLengths,p,opsMovie];
		
		(** plot: intensity **)
		opsMovie={
			AspectRatio->1/2
			,Epilog->{
					Text[Style[p["plot","timeLabel"],p["plot","fontSize"]],Scaled[{0.5,-0.32}]]
					,Rotate[Text[Style[p["plot","intensityLabel"],p["plot","fontSize"]],Scaled[{-0.28,0.5}]],90\[Degree]]
					}
			,ImagePadding->{ipx,{32,8}}
			,ImageSize->55 p["plot","mm"]
			,PlotStyle->AbsoluteThickness[3]
		};
		pI=fig3PlotIntensity[timeData,iData,stepStart,stepEnd,ipx,ops,tTicks,tTicks2,showTumbleShadingQ,tumbleStarts,tumbleLengths,p,opsMovie,1.2];
		
		(** plot: intensity spectrum **)
		pSpectrum=siMoviePlotIntensitySpectrum[t,freqData,trajDataAll,ageString,p,ipx,runColor];
		
		
		(** compose frame **)
		plot=Grid[{
			{timeLabel,SpanFromLeft},
			{mechanicsLabel,opticsLabel},
			{
				Show[plotBacteriumSphere,ImageMargins->{{40,0},{0,0}}]
				,Rasterize[Show[intensityImages0[[t]],ImageSize->296,ImageMargins->{{80,0},{0,0}}]]
			}
			,{
				Column[{pRT,pTheta,pPhi},Alignment->Center,Spacings->0]
				,Column[{pI,pSpectrum},Alignment->Center,Spacings->0]
			}
		},Spacings->0];
		
		(** export frame **)
		Export[FileNameJoin[{outputDirectory,"p_"<>IntegerString[counter++,10,IntegerLength[Length[trajectory]]+1]<>".png"}],plot];
	
	,{t,1,Length[trajectory],frameStep}];
]


(* ::Subsection::Closed:: *)
(*createPlotDroplet*)


(* ::Input:: *)
(*plotDroplet=createPlotDroplet[position];*)


createPlotDroplet[position_]:=Module[{},
	
	(** output **)
	ParametricPlot3D[{
			Evaluate[RotationMatrix[{{0,0,1},position}] . er[\[Theta],\[Phi]]]
			,Evaluate[RotationMatrix[{{0,0,1},position}] . er[\[Pi]-\[Theta],\[Phi]]]
		},{\[Theta],0,0.5\[Pi]},{\[Phi],0,2\[Pi]}
		,Axes->False
		,Boxed->False
		,Lighting->"Neutral"
		,Mesh->None
		,PlotPoints->50
		,PlotStyle->Evaluate[Directive[Opacity[0.2],#,Specularity[White,65]]&/@{Lighter@Red,Lighter[Gray,0.7]}]
	]
]


(* ::Subsection::Closed:: *)
(*createPlotSphericalCap*)


createPlotSphericalCap[sphereRadius_,\[Theta]PlotMaxDeg_]:=ParametricPlot3D[sphereRadius er[\[Theta],\[Phi]],{\[Theta],0,\[Theta]PlotMaxDeg \[Pi]/180},{\[Phi],0,2\[Pi]}
	,Axes->False
	,BoundaryStyle->Directive[Black,AbsoluteThickness[0.5 p["plot","lineThickness"]]]
	,Boxed->False
	,Lighting->"Neutral"
	,Mesh->{1,5}
	,MeshStyle->Directive[Black,AbsoluteThickness[0.5 p["plot","lineThickness"]]]
	,PlotPoints->100
	,PlotStyle->{Directive[Opacity[0.9],Lighter[Gray,0.5],Specularity[White,50]]}
]


(* ::Subsection::Closed:: *)
(*fig3*)


(* ::Input:: *)
(*fig3[p];*)


(* ::Subsubsection::Closed:: *)
(*prepare data*)


(* ::Input:: *)
(*(** for schematic Figs. 3a,b **)*)
(*calcGeometricalParameters[];*)
(*{\[Theta]Bac,\[Phi]Bac}={1.0\[Pi]/4,-1.05\[Pi]/4};*)
(*position=(rEcoliShort+rDroplet)er[\[Theta]Bac,\[Phi]Bac];*)
(*bacDirector=rEcoliShort Normalize[1 e\[Theta][\[Theta]Bac,\[Phi]Bac]+1 e\[Phi][\[Theta]Bac,\[Phi]Bac]];*)
(*plotDroplet=createPlotDroplet[position];*)
(**)
(*pthOutFig3=FileNameJoin[{p["path","figs"],"figure 3"}];*)
(*If[!DirectoryQ[#],CreateDirectory[#]]&@pthOutFig3;*)


(* ::Subsubsection::Closed:: *)
(*Fig. 3a*)


(* ::Input:: *)
(*plot=fig3aPlot[plotDroplet,\[Theta]Bac,\[Phi]Bac,position,bacDirector,rEcoliShort]*)
(*Export[FileNameJoin[{pthOutFig3,"fig3a.png"}],plot,Background->None,ImageResolution->300]*)


(* ::Subsubsection::Closed:: *)
(*Fig. 3b*)


(* ::Input:: *)
(*plot=fig3bPlot[plotDroplet,\[Theta]Bac,\[Phi]Bac,position,bacDirector,rEcoliShort]*)
(*Export[FileNameJoin[{pthOutFig3,"fig3b.png"}],plot,Background->None,ImageResolution->300]*)


(* ::Subsubsection::Closed:: *)
(*code*)


fig3[p_]:=Module[{plotDroplet,\[Theta]Bac,\[Phi]Bac,position,bacDirector,fig3a,fig3b},
	
	(** for schematic Figs. 3a,b **)
	{{rDropletDim,rEcoliShortDim,rEcoliLongDim,rInterfaceDim,\[Rho]FDim,\[Rho]HDim,dDim,mDim,rcomDim,volumeRatio},{rDroplet,rEcoliShort,rEcoliLong,rInterface,\[Rho]F,\[Rho]H,d,m,rcm}}=calcGeometricalParameters[];
	{\[Theta]Bac,\[Phi]Bac}={1.0\[Pi]/4,-1.05\[Pi]/4};
	position=(rEcoliShort+rDroplet)er[\[Theta]Bac,\[Phi]Bac];
	bacDirector=rEcoliShort Normalize[1 e\[Theta][\[Theta]Bac,\[Phi]Bac]+1 e\[Phi][\[Theta]Bac,\[Phi]Bac]];
	plotDroplet=createPlotDroplet[position];
	
	pthOutFig3=FileNameJoin[{p["path","figs"],"figure 3"}];
	If[!DirectoryQ[#],CreateDirectory[#]]&@pthOutFig3;
	
	fig3a=fig3aPlot[plotDroplet,\[Theta]Bac,\[Phi]Bac,position,bacDirector,rEcoliShort];
	Export[FileNameJoin[{pthOutFig3,"fig3a.png"}],fig3a,Background->None,ImageResolution->300];
	Print[fig3a];
	
	fig3b=fig3bPlot[plotDroplet,\[Theta]Bac,\[Phi]Bac,position,bacDirector,rEcoliShort];
	Export[FileNameJoin[{pthOutFig3,"fig3b.png"}],fig3b,Background->None,ImageResolution->300];
	Print[fig3b];
	
	
	
]


(* ::Subsection::Closed:: *)
(*fig3aPlot*)


(* ::Input:: *)
(*plot=fig3aPlot[\[Theta]Bac,\[Phi]Bac,position,bacDirector]*)


fig3aPlot[plotDroplet_,\[Theta]Bac_,\[Phi]Bac_,position_,bacDirector_,rEcoliShort_]:=Module[{plotArrows},
	
	plotArrows=Graphics3D[{
		(*{Opacity[0.3],Sphere[]},*)
		(** coordinate frame **)
		{
			Black,
			Arrow[Tube[{{0,0,0},{1.3,0,0}}]],
			Arrow[Tube[{{0,0,0},{0,-1.3,0}}]],
			Arrow[Tube[{{0,0,0},{0,0,1.3}}]]
		},
		(** force balance arrows: cell propulsion force components **)
		Darker[Green,0.6],
		Arrow[Tube[{#,#+0.61e\[Theta][\[Theta]Bac,\[Phi]Bac]}&@position]],
		Arrow[Tube[{#,#+0.61e\[Phi][\[Theta]Bac,\[Phi]Bac]}&@position]],
		(** force balance arrows: cell propulsion force **)
		Darker[Green,0.2],
		Arrow[Tube[{#,#+0.61Normalize[e\[Theta][\[Theta]Bac,\[Phi]Bac]+e\[Phi][\[Theta]Bac,\[Phi]Bac]]}&@position]],
		(** gravitational restoring force **)
		Darker@Orange,
		Arrow[Tube[{#,#-0.61e\[Theta][\[Theta]Bac,\[Phi]Bac]}&@position]],
		(** position vector **)
		Gray,
		Arrow[Tube[{{0,0,0},0.95er[\[Theta]Bac,\[Phi]Bac]}]],
		(** projection of position vector on middle surface **)
		Darker@Gray,
		Line[{{0,0,0},{1,1,0}er[\[Theta]Bac,\[Phi]Bac]}],
		Line[{er[\[Theta]Bac,\[Phi]Bac],{1,1,0}er[\[Theta]Bac,\[Phi]Bac]}],
		circle3D[{0,0,0},1,{0,0,1}],
		(** bacterium **)
		,{
			Darker[Green,0.2]
			,CapsuleShape[{#-bacDirector,#+bacDirector},rEcoliShort]&@position
		}
		(*Sphere[(0.2+3.5)/3.5er[\[Theta]Bac,\[Phi]Bac],0.2/3.5]*)
		}
		,Background->None
		,Boxed->False
		,Lighting->"Neutral"
		,ViewPoint->{0.97067,-2.9757,1.28559}
		,ViewProjection->Automatic
		,ViewRange->All
		,ViewVector->Automatic
		,ViewVertical->{0.0697,-0.23,0.97069}
	];
	
	(** output **)
	Show[plotArrows,plotDroplet]
]


(* ::Subsection::Closed:: *)
(*fig3bPlot*)


Clear[fig3bPlot];
fig3bPlot[plotDroplet_,\[Theta]Bac_,\[Phi]Bac_,position_,bacDirector_,rEcoliShort_,showCoordinateAxes_:True]:=Module[{plotArrows},
	
	plotArrows=Graphics3D[{
		(** coordinate frame **)
		If[showCoordinateAxes,
			{
				Black
				(*,Arrow[Tube[{{0,0,0},{1.3,0,0}}]]*)
				,Arrow[Tube[{{0,0,0},{0,-1.3,0}}]]
				,Arrow[Tube[{{0,0,0},{0,0,1.3}}]]
			}
		,
			Nothing
		],
		(** arrows **)
		Darker[Green,0.6],
		Arrow[Tube[{#,#+0.6e\[Theta][\[Theta]Bac,\[Phi]Bac]}&@position]],
		(*Arrow[Tube[{#,#+0.4e\[Phi][\[Theta]Bac,\[Phi]Bac]}&@er[\[Theta]Bac,\[Phi]Bac]]],*)
		Darker@Orange,
		Arrow[Tube[{#,#-0.6e\[Theta][\[Theta]Bac,\[Phi]Bac]}&@position]],
		Arrow[Tube[{#,#+0.4{0,0,-1}}&@(-0.525er[\[Theta]Bac,\[Phi]Bac])]],
		(*Red,
		Arrow[Tube[{#,#+0.35e\[Theta][\[Theta]Bac,\[Phi]Bac]+0.35e\[Phi][\[Theta]Bac,\[Phi]Bac]}&@er[\[Theta]Bac,\[Phi]Bac]]],*)
		Gray,
		Arrow[Tube[{{0,0,0},er[\[Theta]Bac,\[Phi]Bac]}]],
		Arrow[Tube[{{0,0,0},-0.49er[\[Theta]Bac,\[Phi]Bac]}]],
		Darker@Gray,
		Sphere[-(0.5+0.02)er[\[Theta]Bac,\[Phi]Bac],0.04],
		(*Line[{{0,0,0},{1,1,0}er[\[Theta]Bac,\[Phi]Bac]}],*)
		(*Line[{er[\[Theta]Bac,\[Phi]Bac],{1,1,0}er[\[Theta]Bac,\[Phi]Bac]}],*)
		(*circle3D[{0,0,Cos[\[Theta]Bac]},Sin[\[Theta]Bac],{0,0,1}],*)
		(*circle3D[{0,0,0},1,{0,0,1}],*)
		Arrowheads[0.02]
		(*Arrow@circle3D[{0,0,0},0.5,{0,0,1},{0.0,-\[Phi]Bac}+3/2 \[Pi]],
		Arrow@circle3D[{0,0,0},0.5,-{1,1,0}er[\[Theta]Bac,\[Phi]Bac-\[Pi]/2],5.25\[Pi]+{0,\[Theta]Bac}]*)
		(*,circle3D[{0,0,0},1,{1,1,0}er[\[Theta]Bac,\[Phi]Bac-\[Pi]/2],{0,2\[Pi]}+\[Theta]Bac]
		,circle3D[{0,0,0},1,{1,1,0}er[\[Theta]Bac,\[Phi]Bac-\[Pi]],{0,2\[Pi]}+\[Theta]Bac]*)
		(** bacterium **)
		,{
			Darker[Green,0.2]
			,CapsuleShape[{#-bacDirector,#+bacDirector},rEcoliShort]&@position
		}
		}
		,Boxed->False
		,Lighting->"Neutral"
		,Background->None
		,ViewProjection->Automatic
		,ViewRange->All
		,ViewVector->Automatic
		,ViewPoint->{-3.2561,-3.02004,-0.1269}
		,ViewVertical->{-0.02149,-0.02185,0.99953}
	];

	(** output **)
	Show[plotArrows,plotDroplet]
]


(* ::Subsection::Closed:: *)
(*fig3CalcData*)


(* ::Input:: *)
(*{timeData,freqData,trajDataAll}=fig3CalcData[trajectories,dt,stepsTransient];*)


fig3CalcData[trajectories_,dt_,stepsTransient_]:=Module[
{windowingQ,nPadding,padMemory,timeData,freqData,trajDataAll,xData,yData,zData,pxData,pyData,pzData,sData,
thetaData,phiData,thetaDataWrapped,phiDataWrapped,iData,meanTheta,freqs,psd,psdi},

	(** fourier settings **)
	windowingQ=False;
	nPadding=2Length[trajectories[[1]]];
	padMemory=ConstantArray[0.0,nPadding];

	timeData=Range[1,Length[trajectories[[1,stepsTransient;;]]-1]]dt;

	Print["Runtime: ",AbsoluteTiming[
	(** analyze trajectories after transient **)
	trajDataAll=Table[
		{xData,yData,zData,pxData,pyData,pzData,sData}=trajectory\[Transpose];
		thetaData=ArcCos/@zData;   (** 0 - top pole (z=+1), \[Pi] - bottom pole (z=-1) **)
		phiData=phaseUnwrap[ArcTan[#[[1]],#[[2]]]&/@trajectory[[All,{1,2}]]];
		thetaDataWrapped=-Abs[thetaData-\[Pi]/2.0]+\[Pi]/2.0;
		phiDataWrapped=-Abs[Mod[phiData,\[Pi]/2.0]-\[Pi]/4.0]+\[Pi]/4.0;
		iData=MapThread[thetaPhiIntensityInterpolation[#1,#2]&,{thetaDataWrapped,phiDataWrapped}];
		meanTheta=Mean[thetaData];
		{freqs,psdi}=calculateTemporalFreqSpectrum[dt,iData,windowingQ,nPadding,padMemory][[{1,2}]];
		
		(** output **)
		{psdi,meanTheta,sData,thetaData,phiData,iData,xData,yData,zData,pxData,pyData,pzData}  
	,{trajectory,trajectories[[All,stepsTransient;;]]}];
	][[1]],"\[ThinSpace]s"];
	freqData=freqs;
	
	(** output **)
	{timeData,freqData,trajDataAll}
]


(* ::Subsection::Closed:: *)
(*fig3CommonPlotOptions*)


fig3CommonPlotOptions[stepPlot_,sData_,dt_,runColor_,tumbleColor_,p_]:=Module[{ipx,colorFunctionRT,trajectoryColoring,ops},
	ipx={54,12};
	colorFunctionRT=getColorFunctionRT[stepPlot,sData,dt,runColor,tumbleColor];
	trajectoryColoring=If[colorizeTumblesQ,
		{
			ColorFunctionScaling->False
			,ColorFunction->colorFunctionRT
			,PlotStyle->AbsoluteThickness[p["plot","lineThickness"]]
		}
	,
		PlotStyle->Directive[AbsoluteThickness[p["plot","lineThickness"]],runColor]
	];
	ops={
		Evaluate[p["plot","ops"]]
		,trajectoryColoring
		,AspectRatio->1/2
		,ImagePadding->{ipx,{30,6}}
		,PlotRangePadding->{Scaled[0.01],Automatic}
	};
	
	(** output **)
	{ipx,ops}
]


(* ::Subsection::Closed:: *)
(*fig3IntensityImageInsets*)


(* ::Input:: *)
(*timeDuration=20;*)
(*{tStart,tEnd}={0,timeDuration};*)
(*{stepStart,stepEnd}={tStart,tEnd}/dt+1;*)
(**)
(*intensityImagesRow=fig3IntensityImageInsets[thetaData,phiData,stepStart,stepEnd,p]*)


(* ::Subsubsection:: *)
(*debug: get intensity images*)


(* ::Input:: *)
(*(** get intensity images **)*)
(*nIntensityImages=4;*)
(*intensityImages0=getIntensityImages[thetaData,phiData,{stepStart,stepEnd},p,nIntensityImages];*)
(*intensityImages=Show[#,Graphics[{Red,Disk[{600,600},80]}],ImagePadding->{{0,0},{0,0}},PlotRangePadding->None,ImageSize->{Automatic,45}]&/@intensityImages0*)


(* ::Subsubsection:: *)
(*debug: visualization test with 3d polygon planes*)


(* ::Input:: *)
(*(** get intensity images **)*)
(*nIntensityImages=6;*)
(*intensityImages0=getIntensityImages[thetaData,phiData,{stepStart,stepEnd},p,nIntensityImages];*)
(**)
(*Graphics3D[{*)
(*Flatten@MapThread[*)
(*{Texture[#1],Polygon[{{#2,-0.5,-0.5},{#2,0.5,-0.5},{#2,0.5,0.5},{#2,-0.5,0.5}},VertexTextureCoordinates->{{0,0},{1,0},{1,1},{0,1}}]}&*)
(*,{intensityImages0,Subdivide[0,2,nIntensityImages-1]}]*)
(*}*)
(*,Boxed->False*)
(*,Lighting->"Neutral"]*)


(* ::Subsubsection:: *)
(*code*)


fig3IntensityImageInsets[thetaData_,phiData_,stepStart_,stepEnd_,nIntensityImages_,p_]:=Module[{intensityImages0,showMarkersQ,intensityImages,markersColors},
	
	(** get intensity images **)
	intensityImages0=getIntensityImages[thetaData,phiData,{stepStart,stepEnd},p,nIntensityImages];
	
	(** resize images **)
	showMarkersQ=True;
	If[showMarkersQ,
		markersColors=timeColormap/@Subdivide[0,1,nIntensityImages-1];
		intensityImages=MapThread[Show[#1,Graphics[{#2,Disk[{600,600},80]}],ImagePadding->{{0,0},{0,0}},PlotRangePadding->None,ImageSize->{Automatic,45}]&,{intensityImages0,markersColors}];
	,
		intensityImages=Show[#,ImagePadding->{{0,0},{0,0}},PlotRangePadding->None,ImageSize->{Automatic,45}]&/@intensityImages0
	];
	
	(** output **)
	Grid[{#&/@intensityImages},Spacings->{0.1,0}]
]



(* ::Subsection::Closed:: *)
(*fig3PlotIntensity*)


(* ::Input:: *)
(*hexToRGB/@{"#ED1C24","#FBB040","#27AAE1","#A96DAE"}*)


(* ::Input:: *)
(*pI=fig3PlotIntensity[timeData,iData,stepStart,stepEnd,ipx,ops,tTicks,tTicks2,showTumbleShadingQ,tumbleStarts,tumbleLengths,p];*)


Clear[fig3PlotIntensity];
fig3PlotIntensity[timeData_,iData_,stepStart_,stepEnd_,ipx_,ops_,tTicks_,tTicks2_,showTumbleShadingQ_,tumbleStarts_,tumbleLengths_,p_,opsMovie_:{},iMaxFactor_:2]:=Module[
{intensityData0,intensityData,iTicks,iTicks2,imin,imax,plotMaxI,markersColors,markersSteps,nMarkers},
	intensityData0={timeData,iData}\[Transpose][[stepStart;;stepEnd]];
	intensityData={#[[1]],Rescale[#[[2]],{imin,imax}]}\[Transpose]&@(intensityData0\[Transpose]);
	iTicks=Table[i,{i,0,1,1}];
	iTicks2={#,""}&/@iTicks;
	{imin,imax}=MinMax[Table[thetaPhiIntensityInterpolation[\[Theta],\[Phi]],{\[Theta],Subdivide[0,\[Pi]/2,100-1]},{\[Phi],Subdivide[0,\[Pi]/4,100-1]}]](*{500,2500}*);     (** imax needs to be large enough to leave space for intensity images in top plot region **)
	plotMaxI=imax*iMaxFactor;  (** boost max plot range by factor 2 so that insets fit **)
	
	nMarkers=4;
	markersColors=timeColormap/@Subdivide[0,1,nMarkers-1];
	markersSteps=Round[Subdivide[0.1,0.9,nMarkers-1](stepEnd-stepStart+1)];
	
	(** output **)
	ListLinePlot[intensityData
		,Evaluate[opsMovie]
		,ImagePadding->{ipx,{35,8}}
		,Evaluate[ops]
		,FrameTicks->{None,{tTicks,tTicks2}}
		,PlotRange->{0,plotMaxI}
		,PlotRangeClipping->False
		,If[showTumbleShadingQ,Prolog->{LightGray,MapThread[Rectangle[{#1,0},{#1+#2,plotMaxI}]&,{tumbleStarts,tumbleLengths}]},{Nothing}]
		,Epilog->{
			(** frame labels **)
			Text[Style[p["plot","timeLabel"],p["plot","fontSize"]],Scaled[{0.5,-0.45}]]
			,Rotate[Text[Style[p["plot","intensityLabel"],p["plot","fontSize"]],Scaled[{-0.17,0.5}]],90\[Degree]]
			(** color markers **)
			,{AbsolutePointSize[4],Table[{markersColors[[i]],Point[intensityData[[markersSteps[[i]]]]]},{i,Length[markersColors]}]}
		}
	]
]


(* ::Subsection::Closed:: *)
(*fig3PlotIntensitySpectrum*)


(* ::Text:: *)
(*plot power spectral density of light intensity (PSD)*)
(*amplitude rescaled according to fit maximum*)


Clear[fig3PlotIntensitySpectrum];
fig3PlotIntensitySpectrum[freqData_,trajDataAll_,ageString_,p_,ipx_,runColor_,opsMovie_:{},saveFig4DataQ_:True]:=Module[
{mPSDi,fMaxIndex,indexStart,indexEnd,iFit,iAmplitude,iWidth,fwhmPoint,mPSDiData,filename,fTicks,fTicks2,pTicks,pTicks2,plot},
	
	(** fit mean power spectral density **)
	mPSDi=Mean[trajDataAll[[All,1]]];
	fMaxIndex=FirstPosition[freqData,_?(#>=p["fourier","fMax"]&)][[1]];
	{indexStart,indexEnd}={1,fMaxIndex};
	iFit=NonlinearModelFit[{freqData,#/Max[#]}\[Transpose][[indexStart;;indexEnd]],{aa/(1+(\[Omega]/bb)^2),aa>0,bb>0},{{aa,1},{bb,0.1}},\[Omega],MaxIterations->1000]&@mPSDi;   (** Lorentz with peak at 0 Hz **)
	{iAmplitude,iWidth}=iFit["BestFitParameters"][[All,2]];
	fwhmPoint={iWidth,0.5};
	mPSDiData={freqData,#/(Max[#]iAmplitude)}\[Transpose][[2;;fMaxIndex]]&@mPSDi;

	(** save data for Fig. 4 **)
	If[saveFig4DataQ,
		filename="meanIntensitySpectrum"<>Capitalize[ageString]<>".hdf5";
		Export[FileNameJoin[{p["path","intensitySpectraNumerical"],filename}],mPSDiData,"HDF5"];
	];

	(** plot: intensity spectrum (PSD)**)
	fTicks=Table[If[Mod[f,0.5]==0,{f,DecimalForm[f,{3,1}],{0.02,0}},{f,"",{0.01,0}}],{f,0,1,0.1}];
	fTicks2={#[[1]],"",#[[3]]}&/@fTicks;
	pTicks=Table[{t,DecimalForm[t,{2,1}]},{t,0,1,0.5}];
	pTicks2={#[[1]],""}&/@pTicks;
	
	plot=Show[
		ListLinePlot[mPSDiData
			,Evaluate[p["plot","ops"]]
			,AspectRatio->1/2
			,PlotRange->{{-0.01,1.01}p["fourier","fMax"],{-0.05,1.25}}
			,FrameTicks->{{pTicks,pTicks2},{fTicks,fTicks2}}
			,ImagePadding->{ipx,{35,8}}
			,PlotStyle->Directive[runColor,AbsoluteThickness[p["plot","lineThickness"]]]
		]
		(*,Plot[iFit[f]/iAmplitude,{f,0,fMax},PlotStyle->Directive[Red,Opacity[0.5]],PlotRange->All]*)
		,PlotRangeClipping->False  (** needed to set frame labels via epilog **)
		,Epilog->{
			(** fwhm highlight **)
			{Red,AbsolutePointSize[5],Point[#],AbsoluteThickness[p["plot","lineThickness"]],Dashing[0.02],HalfLine[{#,{0,0.5}}],HalfLine[{#,{#[[1]],0}}]}&@fwhmPoint,
			(** frame labels **)
			{
				Text[Style[p["plot","frequencyLabel"],p["plot","fontSize"]],Scaled[{0.5,-0.47}]]
				,Rotate[Text[Style[p["plot","spectrumLabel"],p["plot","fontSize"]],Scaled[{-0.37,0.5}]],90\[Degree]]
			}
		}
	];
	
	(** output **)
	{plot,fwhmPoint}
]


(* ::Subsection::Closed:: *)
(*fig3PlotPhi*)


(* ::Input:: *)
(*pPhi=fig3PlotPhi[timeData,phiData,stepStart,stepEnd,\[Phi]OffsetTimeseries,ipx,ops,tTicks,tTicks2,showTumbleShadingQ,tumbleStarts,tumbleLengths,p]*)


Clear[fig3PlotPhi];
fig3PlotPhi[timeData_,phiData_,stepStart_,stepEnd_,\[Phi]OffsetTimeseries_,ipx_,ops_,tTicks_,tTicks2_,showTumbleShadingQ_,tumbleStarts_,tumbleLengths_,p_,opsMovie_:{}]:=Module[{\[Phi]OffsetMean,\[Phi]Ticks,\[Phi]Ticks2},
	\[Phi]OffsetMean=\[Pi]-Mean[phiData[[stepStart;;stepEnd]]];
	\[Phi]Ticks=Table[\[Phi],{\[Phi],0,360,360}];
	\[Phi]Ticks2={#,""}&/@\[Phi]Ticks;
	
	(** output **)
	ListLinePlot[{timeData,180/\[Pi] Mod[phiData+\[Phi]OffsetMean+\[Phi]OffsetTimeseries,2\[Pi]]}\[Transpose][[stepStart;;stepEnd]]
		,Evaluate[opsMovie]
		,ImagePadding->{ipx,{5,8}}
		,AspectRatio->1/4.5
		,Evaluate[ops]
		,FrameTicks->{{\[Phi]Ticks,\[Phi]Ticks2},{tTicks,tTicks2}}
		,PlotRange->{0,360}
		,PlotRangeClipping->False
		,If[showTumbleShadingQ,Prolog->{LightGray,MapThread[Rectangle[{#1,0},{#1+#2,360}]&,{tumbleStarts,tumbleLengths}]},{Nothing}]
		(** frame labels **)
		,Epilog->{
			Rotate[Text[Style["\!\(\*StyleBox[\"\[Phi]\",\nFontSlant->\"Italic\"]\) (\[Degree])",p["plot","fontSize"]],Scaled[{-0.33,0.5}]],90\[Degree]]
		}
	]
]


(* ::Subsection::Closed:: *)
(*fig3PlotRT*)


(* ::Input:: *)
(*showTumbleShadingQ=True;*)
(*pRT=fig3PlotRT[timeData,sData,stepStart,stepEnd,ipx,ops,tTicks2,showTumbleShadingQ,tumbleStarts,tumbleLengths,p]*)


(* ::Subsubsection:: *)
(*code*)


Clear[fig3PlotRT];
fig3PlotRT[timeData_,sData_,stepStart_,stepEnd_,ipx_,ops_,tTicks2_,showTumbleShadingQ_,tumbleStarts_,tumbleLengths_,p_,opsMovie_:{}]:=Module[
{sTicks,sTicks2,tsdata,runStates,tumbleStates},
	
	sTicks={{0,"Tumble"},{1,"Run  "}};
	sTicks2={#[[1]],""}&/@sTicks;
	
	(** -(x-1) flips values 0 and 1, so that Run=1 and Tumble=0 for a nicer plot label; TODO swap meaning in C++ code and elsewhere **)
	tsdata={timeData,-(sData-1)}\[Transpose][[stepStart;;stepEnd]];
	runStates=Select[tsdata,#[[2]]<0.5&];
	tumbleStates=Select[tsdata,#[[2]]>0.5&];
	If[Length[runStates]==0,Print["Warning: No run states!"]; runStates=Nothing;];
	If[Length[tumbleStates]==0,Print["Warning: No tumble states!"]; tumbleStates=Nothing;];
	
	(** output **)
	ListPlot[{runStates,tumbleStates}
		,Evaluate[opsMovie]
		,ImagePadding->{ipx,{5,8}}
		,PlotRangePadding->{Automatic,0.2}
		,Evaluate[ops]
		,FrameTicks->{{sTicks,sTicks2},{tTicks2,tTicks2}}
		,Prolog->If[showTumbleShadingQ,{LightGray,MapThread[Rectangle[{#1,-1},{#1+#2,2}]&,{tumbleStarts,tumbleLengths}]},{Nothing}]
	]
]


(* ::Subsection::Closed:: *)
(*fig3PlotState*)


(* ::Input:: *)
(*ageString="young";(*"old"*)*)
(*trajIndex=1;*)
(*{showTumbleShadingQ,colorizeTumblesQ}={True,True};*)
(*\[Theta]PlotMaxDeg=60;*)
(*\[Phi]OffsetHemisphere=-0.5\[Pi];*)
(*\[Phi]OffsetTimeseries=0.0;*)
(*nIntensityImages=4;*)
(*{fwhmPoint,thetaData,phiData,subplots}=fig3PlotState[ageString,timeData,freqData,trajDataAll,showTumbleShadingQ,colorizeTumblesQ,\[Theta]PlotMaxDeg,dt,pthOut,p,nIntensityImages,trajIndex];*)


(* ::Text:: *)
(*Note: input "trajData" is assumed to already have transient removed*)


(* ::Subsubsection::Closed:: *)
(*debug: RT shading*)


(* ::Input:: *)
(*trajIndex=1;*)
(*{sData,thetaData,phiData,iData}=trajData[[trajIndex,{3,4,5,6}]];*)
(**)
(*timeDuration=20;*)
(*{tStart,tEnd}={0,timeDuration};*)
(*{stepStart,stepEnd}={tStart,tEnd}/dt+1;*)
(**)
(*If[showTumbleShadingQ,*)
(*	rtData=sData[[stepStart;;stepEnd]];*)
(*	tumbleStarts=dt(1+Flatten[Position[Differences@({0}~Join~rtData),1.0]])+tStart;*)
(*	tumbleLengths=dt Length/@Select[Split[rtData],#[[1]]>0.5&];*)
(*];*)
(**)
(*ageString=If[True,"young","old"];*)
(*color=p["color",ageString];*)
(*{runColor,tumbleColor}={color,Black};*)
(*ipx={50,15};*)
(*colorizeTumblesQ=True;*)
(*ops={*)
(*	Evaluate[p["plot","ops"]]*)
(*	,If[colorizeTumblesQ,*)
(*		{*)
(*			ColorFunctionScaling->False*)
(*			,ColorFunction->Function[{t,y},If[sData[[Round[t/dt]]]==0,runColor,tumbleColor]]*)
(*			,PlotStyle->AbsoluteThickness[p["plot","lineThickness"]]*)
(*		}*)
(*	,*)
(*		PlotStyle->Directive[AbsoluteThickness[p["plot","lineThickness"]],runColor]*)
(*	]*)
(*	,AspectRatio->1/2*)
(*	,ImagePadding->{ipx,{30,6}}*)
(*	,PlotRangePadding->{Scaled[0.01],Automatic}*)
(*};*)
(**)
(*ps=ListLinePlot[{timeData,sData}\[Transpose][[stepStart;;stepEnd]]*)
(*	,Evaluate[ops]*)
(*	,FrameLabel->{p["plot","timeLabel"],"RT state"}*)
(*	,If[showTumbleShadingQ,Prolog->{LightGray,MapThread[Rectangle[{#1,0},{#1+#2,1}]&,{tumbleStarts,tumbleLengths}]},Nothing]*)
(*]*)


(* ::Subsubsection::Closed:: *)
(*debug: RT shaded plot quality*)


(* ::Input:: *)
(*showThetaMaxQ=False;*)
(*invertThetaPlotQ=False;*)
(**)
(*tTicks=Table[{t+tStart,t},{t,0,timeDuration,10}];*)
(*tTicks2={#[[1]],""}&/@tTicks;*)
(*\[Theta]Ticks={0,\[Theta]PlotMaxDeg/2,\[Theta]PlotMaxDeg};*)
(*\[Theta]Ticks2={#,""}&/@\[Theta]Ticks;*)
(**)
(**)
(*ListLinePlot[{timeData,180/\[Pi] thetaData}\[Transpose][[stepStart;;stepEnd]]*)
(*	,ImagePadding->{ipx,{10,6}}*)
(*	,Evaluate[ops]*)
(*	,FrameTicks->{{\[Theta]Ticks,\[Theta]Ticks2},{tTicks2,tTicks2}}*)
(*	,ScalingFunctions->If[invertThetaPlotQ,"Reverse",Identity]*)
(*	,PlotRange->{0,\[Theta]PlotMaxDeg}*)
(*	,PlotRangeClipping->False*)
(*	,Prolog->{*)
(*		If[showTumbleShadingQ,{LightGray,MapThread[Rectangle[{#1,0},{#1+#2,\[Theta]PlotMaxDeg}]&,{tumbleStarts,tumbleLengths}]},Nothing]*)
(*		,If[showThetaMaxQ,{Black,Dashed,InfiniteLine[{{0.0,#},{1.0,#}}]&@(180/\[Pi] \[Theta]maxModel)},Nothing]*)
(*	}*)
(*	,Epilog->{Rotate[Text[Style["\!\(\*StyleBox[\"\[Theta]\",\nFontSlant->\"Italic\"]\) (\[Degree])",p["plot","fontSize"]],Scaled[{-0.33,0.5}]],90\[Degree]]}*)
(*]*)


(* ::Subsubsection::Closed:: *)
(*debug: all*)


(* ::Input:: *)
(*youngOldQ=True;*)
(*trajIndex=1;*)
(*{showTumbleShadingQ,colorizeTumblesQ}={True,True};*)
(*\[Theta]PlotMaxDeg=60;*)
(*\[Phi]OffsetHemisphere=-0.5\[Pi];*)
(*\[Phi]OffsetTimeseries=0.0;*)
(*{timeDuration,tStart,stepStart,stepEnd}=fig3SetTime[dt];*)
(*stepPlot=stepEnd-stepStart;*)
(**)
(*{showThetaMaxQ,sData,thetaData,phiData,iData,ageString,runColor,tumbleColor,ipx,trajectoryColoring,ops,tTicks,tTicks2,\[Theta]maxModel,*)
(*rtData,tumbleStarts,tumbleLengths}=fig3PrepareData[stepPlot,trajDataAll,trajIndex,youngOldQ,colorizeTumblesQ,u0Dim,dt,tStart,stepStart,stepEnd,timeDuration,showTumbleShadingQ,p];*)
(*	*)
(*(** plot: trajectory on spherical cap **)*)
(*pTrajectory=fig3TrajectoryOnSphericalCap[trajDataAll[[trajIndex]],\[Theta]PlotMaxDeg,stepStart,stepEnd,colorizeTumblesQ,showThetaMaxQ,\[Theta]maxModel,runColor,tumbleColor,rtData,\[Phi]OffsetHemisphere,rDroplet,rEcoliShort];*)
(*	*)
(*(** plot: RT state **)*)
(*pRT=fig3PlotRT[timeData,sData,stepStart,stepEnd,ipx,ops,tTicks2,showTumbleShadingQ,tumbleStarts,tumbleLengths,p]*)


(* ::Subsubsection:: *)
(*code*)


Clear[fig3PlotState];
fig3PlotState[ageString_,timeData_,freqData_,trajDataAll_,showTumbleShadingQ_,colorizeTumblesQ_,\[Theta]PlotMaxDeg_,dt_,pthOut_,p_,nIntensityImages_:4,trajIndex_:1,\[Phi]OffsetHemisphere_:0.0,\[Phi]OffsetTimeseries_:0.0]:=Module[
{showThetaMaxQ,sData,thetaData,phiData,iData,fwhmPoint,stepPlot,
runColor,tumbleColor,ops,rtData,tumbleStarts,tumbleLengths,\[Theta]maxModel,ipx,timeDuration,tStart,stepStart,stepEnd,
tTicks,tTicks2,pTrajectory,pRT,pTheta,pPhi,intensityImagesRow,intensityData0,intensityData,pI,pSpectrum,subplots0,subplots,plot},
	
	{timeDuration,tStart,stepStart,stepEnd}=fig3SetTime[dt];
	stepPlot=stepEnd-stepStart;
	
	{showThetaMaxQ,sData,thetaData,phiData,iData,runColor,tumbleColor,ipx,ops,tTicks,tTicks2,\[Theta]maxModel,
	rtData,tumbleStarts,tumbleLengths}=fig3PrepareData[stepPlot,trajDataAll,trajIndex,ageString,colorizeTumblesQ,u0Dim,dt,tStart,stepStart,stepEnd,timeDuration,showTumbleShadingQ,p];
	
	(** plot: trajectory on spherical cap **)
	pTrajectory=fig3TrajectoryOnSphericalCap[trajDataAll[[trajIndex]],\[Theta]PlotMaxDeg,stepStart,stepEnd,colorizeTumblesQ,showThetaMaxQ,\[Theta]maxModel,runColor,tumbleColor,rtData,\[Phi]OffsetHemisphere,rDroplet,rEcoliShort];
	
	(** plot: RT state **)
	pRT=fig3PlotRT[timeData,sData,stepStart,stepEnd,ipx,ops,tTicks2,showTumbleShadingQ,tumbleStarts,tumbleLengths,p];
	
	(** plots: polar and azimuthal angle **)
	pTheta=fig3PlotTheta[\[Theta]PlotMaxDeg,timeData,thetaData,stepStart,stepEnd,ipx,ops,tTicks2,showTumbleShadingQ,tumbleStarts,tumbleLengths,showThetaMaxQ,\[Theta]maxModel,p];
	pPhi=fig3PlotPhi[timeData,phiData,stepStart,stepEnd,\[Phi]OffsetTimeseries,ipx,ops,tTicks,tTicks2,showTumbleShadingQ,tumbleStarts,tumbleLengths,p];
	
	(** plot: intensity images **)
	intensityImagesRow=fig3IntensityImageInsets[thetaData,phiData,stepStart,stepEnd,nIntensityImages,p];
	
	(** plot: intensity **)
	pI=fig3PlotIntensity[timeData,iData,stepStart,stepEnd,ipx,ops,tTicks,tTicks2,showTumbleShadingQ,tumbleStarts,tumbleLengths,p];
	
	(** plot: intensity spectrum **)
	{pSpectrum,fwhmPoint}=fig3PlotIntensitySpectrum[freqData,trajDataAll,ageString,p,ipx,runColor];
	
	(** organize single subplots as an array for easier export **)
	subplots0=Show[#,ImageSize->45 p["plot","mm"],ImageMargins->0]&/@{pTrajectory,pRT,pTheta,pPhi,pI,pSpectrum};
	subplots=subplots0[[;;4]]~Join~{intensityImagesRow}~Join~subplots0[[5;;]];
	plot=Grid[{#}&/@subplots];
	
	(** export composite plot **)
	Export[FileNameJoin[{pthOut,"mma_figure3_"<>ageString<>".png"}],plot,ImageResolution->300];
	
	(** subplots for composing in Adobe Illustrator **)
	(** png export **)
	Do[Export[FileNameJoin[{pthOut,"figure3"<>"_"<>ageString<>"_part_"<>ToString[i]<>".png"}],subplots[[i]],ImageResolution->300],{i,Length[subplots]}];
	(** pdf export **)
	Do[Export[FileNameJoin[{pthOut,"figure3"<>"_"<>ageString<>"_part_"<>ToString[i]<>".pdf"}],subplots[[i]]],{i,Length[subplots]}];
	
	(** show plot **)
	Print[Magnify[plot,1.5]];
	
	(** output **)
	{fwhmPoint,thetaData,phiData,subplots}
]


(* ::Subsection::Closed:: *)
(*fig3PlotTheta*)


(* ::Input:: *)
(*pTheta=fig3PlotTheta[\[Theta]PlotMaxDeg,timeData,thetaData,stepStart,stepEnd,ipx,ops,tTicks2,showTumbleShadingQ,tumbleStarts,tumbleLengths,showThetaMaxQ,\[Theta]maxModel,p]*)


Clear[fig3PlotTheta];
fig3PlotTheta[\[Theta]PlotMaxDeg_,timeData_,thetaData_,stepStart_,stepEnd_,ipx_,ops_,tTicks2_,showTumbleShadingQ_,tumbleStarts_,tumbleLengths_,showThetaMaxQ_,\[Theta]maxModel_,p_,opsMovie_:{}]:=Module[
{invertThetaPlotQ,\[Theta]Ticks,\[Theta]Ticks2},
	invertThetaPlotQ=False;
	\[Theta]Ticks={0,\[Theta]PlotMaxDeg/2,\[Theta]PlotMaxDeg};
	\[Theta]Ticks2={#,""}&/@\[Theta]Ticks;

	(** output **)
	ListLinePlot[{timeData,180/\[Pi] thetaData}\[Transpose][[stepStart;;stepEnd]]
		,Evaluate[opsMovie]
		,ImagePadding->{ipx,{5,8}}
		,Evaluate[ops]
		,FrameTicks->{{\[Theta]Ticks,\[Theta]Ticks2},{tTicks2,tTicks2}}
		,ScalingFunctions->If[invertThetaPlotQ,"Reverse",Identity]
		,PlotRange->{0,\[Theta]PlotMaxDeg}
		,PlotRangeClipping->False
		,Prolog->{
			If[showTumbleShadingQ,{LightGray,MapThread[Rectangle[{#1,0},{#1+#2,\[Theta]PlotMaxDeg}]&,{tumbleStarts,tumbleLengths}]},{Nothing}]
			,If[showThetaMaxQ,{Black,Dashed,InfiniteLine[{{0.0,#},{1.0,#}}]&@(180/\[Pi] \[Theta]maxModel)},{Nothing}]
		}
		,Epilog->{Rotate[Text[Style["\!\(\*StyleBox[\"\[Theta]\",\nFontSlant->\"Italic\"]\) (\[Degree])",p["plot","fontSize"]],Scaled[{-0.33,0.5}]],90\[Degree]]}
	]

]


(* ::Subsection::Closed:: *)
(*fig3PrepareData*)


(* ::Input:: *)
(*{showThetaMaxQ,sData,thetaData,phiData,iData,runColor,tumbleColor,ipx,ops,tTicks,tTicks2,\[Theta]maxModel,rtData,tumbleStarts,tumbleLengths}=fig3PrepareData[stepPlot,trajDataAll,trajIndex,youngOldQ,colorizeTumblesQ,u0Dim,dt,tStart,stepStart,stepEnd,timeDuration,showTumbleShadingQ,p];*)


(* ::Subsubsection:: *)
(*code*)


fig3PrepareData[stepPlot_,trajDataAll_,trajIndex_,ageString_,colorizeTumblesQ_,u0Dim_,dt_,tStart_,stepStart_,stepEnd_,timeDuration_,showTumbleShadingQ_,p_]:=Module[
{showThetaMaxQ,sData,thetaData,phiData,iData,runColor,tumbleColor,ipx,ops,tTicks,tTicks2,\[Theta]maxModel,
rtData,tumbleStarts,tumbleLengths},
	
	(** options **)
	showThetaMaxQ=False;
	
	(** get state data **)
	{sData,thetaData,phiData,iData}=trajDataAll[[trajIndex,{3,4,5,6}]];
	
	(** options **)
	{runColor,tumbleColor}={#,(*Black*)Darker[#,0.8]}&@p["color",ageString];
	{ipx,ops}=fig3CommonPlotOptions[stepPlot,sData,dt,runColor,tumbleColor,p];
	tTicks=Table[{t+tStart,t},{t,0,timeDuration,10}];
	tTicks2={#[[1]],""}&/@tTicks;
	
	If[showThetaMaxQ,
		\[Theta]maxModel=maxTheta[u0Dim];
		If[Im[\[Theta]maxModel]!=0,Print["Warning: \[Theta]maxModel is complex!"]];
	];
	
	If[showTumbleShadingQ,
		rtData=sData[[stepStart;;stepEnd]];
		tumbleStarts=dt(Flatten[Position[Differences@({0}~Join~rtData),1.0]])+tStart;
		tumbleLengths=dt Length/@Select[Split[rtData],#[[1]]>0.5&];
	];
	
	(** output **)
	{showThetaMaxQ,sData,thetaData,phiData,iData,runColor,tumbleColor,ipx,ops,tTicks,tTicks2,\[Theta]maxModel,rtData,tumbleStarts,tumbleLengths}
]


(* ::Subsection::Closed:: *)
(*fig3SetTime*)


(* ::Input:: *)
(*{timeDuration,tStart,stepStart,stepEnd}=fig3SetTime[dt]*)


(* ::Text:: *)
(*[timeDuration] = s*)
(*transient is already deducted before this*)


Clear[fig3SetTime];
fig3SetTime[dt_,timeDuration_:20]:=Module[{tStart,tEnd,stepStart,stepEnd},
	{tStart,tEnd}={0,timeDuration};
	{stepStart,stepEnd}=Round[{tStart,tEnd}/dt+1]; (** offset of 1 needed b/c mma arrays start at 1 instead of 0 **)
	
	(** output **)
	{timeDuration,tStart,stepStart,stepEnd}
]


(* ::Subsection::Closed:: *)
(*fig3TrajectoryOnSphericalCap*)


(* ::Text:: *)
(*plot: trajectory on spherical cap*)


(* ::Input:: *)
(*pTrajectory=fig3TrajectoryOnSphericalCap[trajData,\[Theta]PlotMaxDeg,stepStart,stepEnd,colorizeTumblesQ,showThetaMaxQ,\[Theta]maxModel,runColor,tumbleColor,rtData,\[Phi]OffsetHemisphere,rDroplet,rEcoliShort]*)


(* ::Subsubsection:: *)
(*debug*)


(* ::Input:: *)
(*timeDuration=20;*)
(*{tStart,tEnd}={0,timeDuration};*)
(*{stepStart,stepEnd}={tStart,tEnd}/dt+1;*)
(*showThetaMaxQ=False;*)
(*\[Theta]maxModel;*)
(*ageString=If[youngOldQ,"young","old"];*)
(*{runColor,tumbleColor}={p["color",ageString],Black};*)
(*trajData=trajDataAll[[trajIndex]];*)
(*{sData,thetaData,phiData,iData}=trajData[[{3,4,5,6}]];*)
(*rtData=sData[[stepStart;;stepEnd]]*)
(**)
(*fig3SphericalCap[trajData,\[Theta]PlotMaxDeg,stepStart,stepEnd,showThetaMaxQ,\[Theta]maxModel,runColor,tumbleColor,rtData,\[Phi]OffsetHemisphere,rDroplet,rEcoliShort]*)


(* ::Subsubsection:: *)
(*code*)


fig3TrajectoryOnSphericalCap[trajData_,\[Theta]PlotMaxDeg_,stepStart_,stepEnd_,colorizeTumblesQ_,showThetaMaxQ_,\[Theta]maxModel_,runColor_,tumbleColor_,rtData_,\[Phi]OffsetHemisphere_,rDroplet_,rEcoliShort_]:=Module[
{sphereRadius,plotCap,plotRangeRadiusMax,rotaMatrix,rotatedTrajectory,bacteriumPosition,bacteriumDirector},
	
	
	sphereRadius=0.99;
	plotCap=createPlotSphericalCap[sphereRadius,\[Theta]PlotMaxDeg];
	plotRangeRadiusMax=1.1 Sin[\[Theta]PlotMaxDeg Degree];
	rotaMatrix=RotationMatrix[\[Phi]OffsetHemisphere,{0,0,1}];
	rotatedTrajectory=1.01 sphereRadius (rotaMatrix . #)&/@Transpose[trajData[[{7,8,9},stepStart;;stepEnd]]];
	bacteriumPosition=(rDroplet+rEcoliShort)rotatedTrajectory[[-1]];
	bacteriumDirector=rEcoliShort rotaMatrix . trajData[[{10,11,12},stepEnd]];

	(** output **)
	Show[
		plotCap,
		Graphics3D[{
			(** mark max angle from model **)
			If[showThetaMaxQ,
				{
					Gray
					,Dashed
					,circle3D[{0,0,sphereRadius Cos[#]},sphereRadius Sin[#]]&@\[Theta]maxModel
				}
				,Nothing
			]
			(** trajectory **)
			,{
				AbsoluteThickness[2p["plot","lineThickness"]]
				,runColor
				,Line[rotatedTrajectory,If[colorizeTumblesQ,VertexColors->Evaluate[If[#==0,runColor,tumbleColor]&/@rtData],Nothing]
				]
			}
			(** bacterium **)
			,{
				Opacity[0.4],runColor,CapsuleShape[{#-bacteriumDirector,#+bacteriumDirector},rEcoliShort]&@bacteriumPosition
			}
			(** z-axis **)
			,{
				Black,Arrowheads[0.03],Arrow[Tube[{{0,0,sphereRadius},{0,0,1.23}},0.003]]
			}
		}]
		,Axes->False
		,Boxed->False
		,Method->{"ShrinkWrap"->True}(*,ImagePadding->{{50,15},{10,10}}*)
		,PlotRange->{plotRangeRadiusMax{-1,1},plotRangeRadiusMax{-1,1},{0.2,1.3}}
		,ViewPoint->{0.398999,-2.9706995,1.57011}
		,ViewVertical->{0.03589,-0.2831,0.9584}
	]
]


(* ::Subsection::Closed:: *)
(*getColorFunctionRT*)


(* ::Input:: *)
(*colorFunctionRT=getColorFunctionRT[step,sData,dt,runColor,tumbleColor];*)


getColorFunctionRT[step_,sData_,dt_,runColor_,tumbleColor_]:=Function[{t,y},If[Round[t/dt]<=step,If[sData[[Round[t/dt]]]==0,runColor,tumbleColor],LightGray]]


(* ::Subsection:: *)
(*getSimulationData*)


(* ::Input:: *)
(*{timeData,freqData,trajDataAll}=getSimulationData[True,u0Dim,runTime0,tumbleTime0,frictionFactor,loadDataQ,p]*)


Clear[getSimulationData];
getSimulationData[youngOldQ_,u0Dim_,runTime0_,tumbleTime0_,frictionFactor_,loadDataQ_,p_]:=Module[{debugPlotQ,plotNoiseQ,detailedForceAnalysisQ,diffRotDim,rDropletDim,
diffS,timeFactor,asymmetry,tfinalSim,nParallelSystems,noiseSeed,tTransient,hetQ,pathLookupMap,swimMultiplier,stepSize,fileName,
rtRate,trRate,swimFactor,gravFactor,args,timeData,freqData,trajDataAll},
	{debugPlotQ,plotNoiseQ,detailedForceAnalysisQ}={False,False,False};
	
	(** parameter values **)
	{diffRotDim,rDropletDim,swimMultiplier}={3.5,3.5*^-6,1.0};
	{diffS,timeFactor,asymmetry}={115,0.186,16.0}; 
	{tfinalSim,nParallelSystems,noiseSeed}={200,200,0};
	tTransient=Max[1,Round[0.2tfinalSim]];
	hetQ=0;
	stepSize=1;
	
	pathLookupMap=FileNameJoin[{p["path","dropletOpticsNumerical"],"lookup_map"}];
	thetaPhiIntensityInterpolation=updateIntensityLookupMap[pathLookupMap];
	
	{rtRate,trRate}=1.0/{runTime0,tumbleTime0};
	calcPhysicalParameters[u0Dim,rDropletDim];
	{swimFactor,gravFactor}=calcForceFactors[rDropletDim,u0Dim,swimMultiplier];
	
	If[!loadDataQ,
		(** run simulation **)
		Print["Running simulation"];
		args={u0Dim,diffRotDim,diffS,timeFactor,rDropletDim,asymmetry,swimMultiplier,tfinalSim,nParallelSystems,noiseSeed,hetQ,rtRate,trRate,frictionFactor};
		runSimulation[args];
		{nSteps,dt,times,trajectories,stepsTransient}=readSimulationOutput[tfinalSim,tTransient,stepSize];
	,
		(** load previous simulation file **)
		fileName=FileNameJoin[{p["path","runTumbleTrajectoriesNumerical"],"runTumbleOutput_"<>If[youngOldQ,"young","old"]<>".bin"}];
		{nSteps,dt,times,trajectories,stepsTransient}=readSimulationOutput[tfinalSim,tTransient,stepSize,fileName,7,10,0.001];
	];

	{timeData,freqData,trajDataAll}=fig3CalcData[trajectories,dt,stepsTransient];
	
	(** output **)
	{timeData,freqData,trajDataAll}
]


(* ::Subsection::Closed:: *)
(*saveNewSimulationResultAsDefault*)


(* ::Subsubsection:: *)
(*code*)


(** Note: only run this when you are sure you want to keep the results! **)
saveNewSimulationResultAsDefault[youngOldString_,p_]:=Module[{oldPthFn,newPthFn},
	If[!MemberQ[{"young","old"},youngOldString],Print["Error: youngOldString ("<>youngOldString<>") must be either \"young\" or \"old\"!"];];
	oldPthFn=FileNameJoin[{NotebookDirectory[],"build","runTumbleOutput.bin"}];
	newPthFn=FileNameJoin[{p["path","runTumbleTrajectoriesNumerical"],"runTumbleOutput_"<>youngOldString<>".bin"}];
	If[!DirectoryQ@#,CreateDirectory@#]&@p["path","runTumbleTrajectoriesNumerical"];
	If[FileExistsQ@#,DeleteFile[#]]&@newPthFn;
	CopyFile[oldPthFn,newPthFn];
]


(* ::Subsection::Closed:: *)
(*siMoviePlotIntensitySpectrum*)


Clear[siMoviePlotIntensitySpectrum];
siMoviePlotIntensitySpectrum[step_,freqData_,trajDataAll_,ageString_,p_,ipx_,runColor_]:=Module[
{mPSDi,fMaxIndex,indexStart,indexEnd,iFit,iAmplitude,iWidth,fwhmPoint,mPSDiData,filename,fTicks,fTicks2,pTicks,pTicks2,plot},
	
	(** perform PSD **)
	
	(** fit mean power spectral density **)
	mPSDi=Mean[trajDataAll[[All,1]]];
	fMaxIndex=FirstPosition[freqData,_?(#>=p["fourier","fMax"]&)][[1]];
	{indexStart,indexEnd}={1,fMaxIndex};
	iFit=NonlinearModelFit[{freqData,#/Max[#]}\[Transpose][[indexStart;;indexEnd]],{aa/(1+(\[Omega]/bb)^2),aa>0,bb>0},{{aa,1},{bb,0.1}},\[Omega],MaxIterations->1000]&@mPSDi;   (** Lorentz with peak at 0 Hz **)
	{iAmplitude,iWidth}=iFit["BestFitParameters"][[All,2]];
	fwhmPoint={iWidth,0.5};
	mPSDiData={freqData,#/(Max[#]iAmplitude)}\[Transpose][[2;;fMaxIndex]]&@mPSDi;

	(** plot: intensity spectrum (PSD)**)
	fTicks=Table[If[Mod[f,0.5]==0,{f,DecimalForm[f,{3,1}],{0.02,0}},{f,"",{0.01,0}}],{f,0,1,0.1}];
	fTicks2={#[[1]],"",#[[3]]}&/@fTicks;
	pTicks=Table[{t,DecimalForm[t,{2,1}]},{t,0,1,0.5}];
	pTicks2={#[[1]],""}&/@pTicks;
	
	plot=Show[
		ListLinePlot[mPSDiData
			,AspectRatio->1/2
			,ImagePadding->{ipx,{36,8}}
			,ImageSize->55p["plot","mm"]
			,Evaluate[p["plot","ops"]]
			,PlotRange->{{-0.01,1.01}p["fourier","fMax"],{-0.05,1.25}}
			,FrameTicks->{{pTicks,pTicks2},{fTicks,fTicks2}}
			,ImagePadding->{ipx,{35,8}}
			,PlotStyle->Directive[runColor,AbsoluteThickness[p["plot","lineThickness"]]]
		]
		(*,Plot[iFit[f]/iAmplitude,{f,0,fMax},PlotStyle->Directive[Red,Opacity[0.5]],PlotRange->All]*)
		,PlotRangeClipping->False
		,Epilog->{
			(** fwhm highlight **)
			(*{Red,AbsolutePointSize[5],Point[#],AbsoluteThickness[p["plot","lineThickness"]],Dashing[0.02],HalfLine[{#,{0,0.5}}],HalfLine[{#,{#[[1]],0}}]}&@fwhmPoint,*)
			(** frame labels **)
			{
				Text[Style[p["plot","frequencyLabel"],p["plot","fontSize"]],Scaled[{0.5,-0.35}]]
				,Rotate[Text[Style["Intensity spectrum",p["plot","fontSize"]],Scaled[{-0.27,0.35}]],90\[Degree]]
			}
		}
	];
	
	(** output **)
	plot
]


(* ::Section:: *)
(*Subfigures*)


(* ::Subsection:: *)
(*Figure 3a,b*)


(* ::Input:: *)
(*fig3[p];*)


(* ::Subsection::Closed:: *)
(*find numerical parameters matching the experiment*)


(* ::Input:: *)
(*(** get experimental width **)*)
(*widthsEx=Flatten[Import[FileNameJoin[{p["path","widthsExperiment"],"experimental_widths.dat"}]]];*)
(**)
(*(** get simulation scan data **)*)
(*dataScan=loadSimulationScanData[p];*)
(**)
(*(** get discretized contours **)*)
(*{widthInterpolationLogLogLin,angleInterpolationLogLogLin}=getScanDataInterpolantsLogLogLin[dataScan];*)
(*{wMin,wMax}={Floor[#[[1]],0.01],Ceiling[#[[2]],0.01]}&@MinMax[widthsEx];*)


(* ::Input:: *)
(*parameterYoung={Log10@#[[2]],Log10@#[[3]],10^6#[[1]]}&@{7.2*^-6,0.9,0.1}*)
(*parameterOld={Log10@#[[2]],Log10@#[[3]],10^6#[[1]]}&@{3.2*^-6,1.6,0.01}*)
(*contourPlot=ContourPlot3D[widthInterpolationLogLogLin[x,y,z],{x,Log10@0.1,Log10@10},{y,Log10@0.01,Log10@1},{z,1,10},Contours->{wMin,wMax},Mesh->None,ContourStyle->Opacity[0.8]];*)
(**)
(*Show[*)
(*contourPlot*)
(*,Graphics3D[{Black,Scale[Sphere[parameterYoung,0.03],{1,1,4}],Scale[Sphere[parameterOld,0.03],{1,1,4}]}]*)
(*]*)


(* ::Subsection:: *)
(*Fig. 3 (left column): young bacterium*)


(* ::Subsubsection::Closed:: *)
(*simulation parameter history*)


(* ::Input:: *)
(*(** 10 um/s - slightly wider; 20 um/s - shorter **)*)
(*{u0Dim,runTime0,tumbleTime0}={10.0*^-6,0.2,0.5};      (** works for young bacterium w=0.18 **)*)
(*{u0Dim,runTime0,tumbleTime0}={10.0*^-6,0.2,0.4};      (** works for young bacterium w=0.19 **)*)
(*{u0Dim,runTime0,tumbleTime0}={10.0*^-6,0.3,0.5};      (** works for young bacterium w=0.17 **)*)
(*{u0Dim,runTime0,tumbleTime0}={10.0*^-6,0.4,0.5};      (** works for young bacterium w=0.14 **)*)
(*{u0Dim,runTime0,tumbleTime0}={10.0*^-6,0.4,0.4};      (** works for young bacterium w=0.12 **)*)
(*{u0Dim,runTime0,tumbleTime0}={10.0*^-6,0.35,0.5};    (** works for young bacterium w=0.15 **)*)
(*{u0Dim,runTime0,tumbleTime0}={10.0*^-6,0.5,0.3};      (** works for young bacterium w=0.11 **)*)
(*{u0Dim,runTime0,tumbleTime0}={10.0*^-6,0.5,0.2};      (** works for young bacterium w=0.10 **)*)
(*{u0Dim,runTime0,tumbleTime0}={10.0*^-6,0.5,0.5};      (** works for young bacterium w=0.12 **)*)
(*{u0Dim,runTime0,tumbleTime0}={10.0*^-6,0.5,0.7};      (** works for young bacterium w=0.13 **)*)
(*{u0Dim,runTime0,tumbleTime0}={10.0*^-6,0.5,0.9};      (** works for young bacterium w=0.15 **)*)
(*{u0Dim,runTime0,tumbleTime0}={8.0*^-6,0.2,0.4};        (** works for young bacterium w=0.17 **)*)
(*{u0Dim,runTime0,tumbleTime0}={7.0*^-6,0.2,0.4};        (** works for young bacterium w=0.16 **)*)
(*{u0Dim,runTime0,tumbleTime0}={6.0*^-6,0.2,0.4};        (** works for young bacterium w=0.15 **)*)
(*{u0Dim,runTime0,tumbleTime0}={10.0*^-6,0.2,0.4};      (** works for young bacterium w=0.23 **)*)
(*{u0Dim,runTime0,tumbleTime0}={10.0*^-6,0.42,0.13};  (** great, but width too small **)*)
(*{u0Dim,runTime0,tumbleTime0}={10.0*^-6,0.1,0.09};    (** found by fit **)*)
(*{u0Dim,runTime0,tumbleTime0}={10.0*^-6,1,0.09};        (** max width **)*)
(*{u0Dim,runTime0,tumbleTime0}={12.0*^-6,0.9,0.1};      (** for friction=0.25 friction0, w=0.29 **)*)
(*{u0Dim,runTime0,tumbleTime0}={10.0*^-6,0.9,0.1};      (** for friction=0.25 friction0, w=0.21 **)*)
(*{u0Dim,runTime0,tumbleTime0}={11.0*^-6,0.9,0.1};      (** for friction=0.25 friction0, w=0.25 **)*)
(*{u0Dim,runTime0,tumbleTime0}={11.0*^-6,0.9,0.1};      (** for friction=0.195 friction0, w=0.27 **)*)
(**)
(*{u0Dim,runTime0,tumbleTime0}={10.5*^-6,0.9,0.1};      (** for friction=0.195 friction0, w=0.25 **)*)
(*{u0Dim,runTime0,tumbleTime0}={10.0*^-6,0.9,0.1};      (** for friction=0.195 friction0, w=0.19 **)*)
(*{u0Dim,runTime0,tumbleTime0}={11.0*^-6,0.9,0.1};*)
(*{u0Dim,runTime0,tumbleTime0}={6.8*^-6,1.0,0.1};*)
(**)


(* ::Subsubsection:: *)
(*run simulation*)


(* ::Input:: *)
(*{u0Dim,runTime0,tumbleTime0,frictionFactor}={7.36*^-6,0.9,0.1,1.0};      (** for friction=0.195 friction0, w=0.23 **)*)
(*loadDataQ=False;*)
(*{timeDataYoung,freqDataYoung,trajDataAllYoung}=getSimulationData[True,u0Dim,runTime0,tumbleTime0,frictionFactor,loadDataQ,p];*)


(* ::Subsubsection:: *)
(*save simulation as new default*)


(* ::Input:: *)
(*(** Note: only run this when you are sure you want to keep the results! **)*)
(*saveNewSimulationResultAsDefault["young",p];*)


(* ::Subsubsection:: *)
(*plot data*)


(* ::Input:: *)
(*nIntensityImages=4;*)
(*trajIndex=1;*)
(*{showTumbleShadingQ,colorizeTumblesQ}={True,True};*)
(*\[Theta]PlotMaxDeg=60;*)
(*\[Phi]OffsetHemisphere=-0.5\[Pi];*)
(*\[Phi]OffsetTimeseries=0.0;*)
(**)
(*{fwhmPointYoung,thetaData,phiData,subplots}=*)
(*fig3PlotState["young",timeDataYoung,freqDataYoung,trajDataAllYoung,showTumbleShadingQ,colorizeTumblesQ,\[Theta]PlotMaxDeg,dt,pthOutFig3,p,nIntensityImages,trajIndex,\[Phi]OffsetHemisphere,\[Phi]OffsetTimeseries];*)
(*fwhmPointYoung*)


(* ::Subsubsection::Closed:: *)
(*test: identify tumbles/runs from trajectory*)


(* ::Input:: *)
(*tumbleSections=Select[{sData,timeData,thetaData}\[Transpose],#[[1]]==1&];*)
(*ListPlot[tumbleSections[[All,3]],ScalingFunctions->"Reverse",Frame->True]*)


(* ::Subsubsection::Closed:: *)
(*SI Video: young bacterium with plots*)


(* ::Input:: *)
(*(** function from droplet optics **)*)
(*solution=solveHemisphereLC[];*)


(* ::Input:: *)
(*(** create video frames **)*)
(*pthFrames=FileNameJoin[{$HomeDirectory,"Desktop","mmamovie_runTumbleSimulation_young_2"}];*)
(*{tmax,tfinalSim}={20,200};*)
(**)
(*Print["movie frame generation took: ",AbsoluteTiming[*)
(*createMovieFramesSimulation[trajectories[[trajIndex,stepsTransient;;stepsTransient+Round[tmax/tfinalSim nSteps]]],"young",dt,trajDataAllYoung,timeDataYoung,freqDataYoung,pthFrames,p,tmax,2*)
(*];*)
(*][[1]],"\[ThinSpace]s"];*)


(* ::Input:: *)
(*(** render video with ffmpeg **)*)
(*fps=50;*)
(*renderVideo[pthFrames,FileNameJoin[{pthVideo,"Supplementary Videos","Supplementary Video 9 Run Tumble Sphere young full.mp4"}],fps];*)


(* ::Subsubsection::Closed:: *)
(*SI Video: only bacterium*)


(* ::Input:: *)
(*(** create video frames **)*)
(*pthFrames=FileNameJoin[{$HomeDirectory,"Desktop","mmamovie_runTumbleSphere_young"}];*)
(*{tmax,tfinalSim}={20,200};*)
(*Print[AbsoluteTiming[*)
(*	createMovieFramesPosition[trajectories[[trajIndex,stepsTransient;;stepsTransient+Round[tmax/tfinalSim nSteps];;2]],1,dt];*)
(*][[1]]];*)
(**)
(*(** render video with ffmpeg **)*)
(*fps=50;*)
(*renderVideo[pthFrames,FileNameJoin[{pthVideo,"Supplementary Videos","Supplementary Video 9 Run Tumble Sphere young.mp4"}],fps];*)


(* ::Input:: *)
(*(** test single frame **)*)
(*createMovieFramesPosition[trajectories[[trajIndex,stepsTransient;;stepsTransient+1;;2]],1,dt]*)


(* ::Subsubsection:: *)
(*SI Video: intensity images*)


(* ::Input:: *)
(*(** get intensity images **)*)
(*nIntensityImages=stepEnd-stepStart+1;*)
(*intensityImages0=getIntensityImages[thetaData,phiData,{stepStart,stepEnd},p,nIntensityImages];*)
(**)
(*(** create video frames **)*)
(*pthFrames=FileNameJoin[{$HomeDirectory,"Desktop","mmamovie_intensity_images_young"}];*)
(*If[!DirectoryQ@#,CreateDirectory@#]&@pthOutput;*)
(*counter=0;*)
(*Do[*)
(*	Export[FileNameJoin[{pthFrames,"p_"<>IntegerString[counter++,10,IntegerLength[Length[intensityImages0]]+1]<>".png"}],intensityImages0[[i]]]*)
(*,{i,Length[intensityImages0]}];*)
(**)
(*(** render video with ffmpeg **)*)
(*fps=50;*)
(*renderVideo[pthFrames,FileNameJoin[{pthVideo,"Supplementary Videos","Supplementary Video 9 intensity young.mp4"}],fps];*)


(* ::Subsubsection::Closed:: *)
(*SI Video: time traces*)


(* ::Input:: *)
(*trajIndex=1;*)
(*{showTumbleShadingQ,colorizeTumblesQ}={True,True};*)
(*\[Theta]PlotMaxDeg=60;*)
(*\[Phi]OffsetHemisphere=-0.5\[Pi];*)
(*\[Phi]OffsetTimeseries=0.0;*)
(*youngOldQ=True;*)
(**)
(*{timeDuration,tStart,stepStart,stepEnd}=fig3SetTime[dt];*)
(*stepPlot=100;*)
(*showTumbleShadingQ=False;*)
(**)
(*{showThetaMaxQ,sData,thetaData,phiData,iData,ageString,runColor,tumbleColor,ipx,ops,tTicks,tTicks2,\[Theta]maxModel,*)
(*rtData,tumbleStarts,tumbleLengths}=fig3PrepareData[stepPlot,trajDataAllYoung,trajIndex,youngOldQ,colorizeTumblesQ,u0Dim,dt,tStart,stepStart,stepEnd,timeDuration,showTumbleShadingQ,p];*)


(* ::Input:: *)
(*nFrames=stepEnd-stepStart;*)
(*counter=0;*)
(*Print["Runtime: ",AbsoluteTiming[*)
(*Do[*)
(*(** update plot option **)*)
(*{ipx,ops}=fig3CommonPlotOptions[stepPlot,sData,dt,runColor,tumbleColor,p];*)
(**)
(*(** plot: RT state **)*)
(*showTumbleShadingQ=False;*)
(*opsMovie={AspectRatio->1/3,PlotStyle->AbsolutePointSize[3]};*)
(*pRT=fig3PlotRT[timeData,sData,stepStart,stepEnd,ipx,ops,tTicks2,showTumbleShadingQ,tumbleStarts,tumbleLengths,p,opsMovie];*)
(**)
(*(** plots: polar and azimuthal angle **)*)
(*opsMovie={AspectRatio->1/3,PlotStyle->AbsoluteThickness[3]};*)
(*pTheta=fig3PlotTheta[\[Theta]PlotMaxDeg,timeData,thetaData,stepStart,stepEnd,ipx,ops,tTicks2,showTumbleShadingQ,tumbleStarts,tumbleLengths,showThetaMaxQ,\[Theta]maxModel,p,opsMovie];opsMovie={AspectRatio->1/3*)
(*,ImagePadding->{ipx,{40,8}}*)
(*,Epilog->{*)
(*		Text[Style[p["plot","timeLabel"],p["plot","fontSize"]],Scaled[{0.5,-0.75}]]*)
(*		,Rotate[Text[Style["\!\(\*StyleBox[\"\[Phi]\",\nFontSlant->\"Italic\"]\) (\[Degree])",p["plot","fontSize"]],Scaled[{-0.33,0.5}]],90\[Degree]]*)
(*}*)
(*,PlotStyle->AbsoluteThickness[3]*)
(*};*)
(*pPhi=fig3PlotPhi[timeData,phiData,stepStart,stepEnd,\[Phi]OffsetTimeseries,ipx,ops,tTicks,tTicks2,showTumbleShadingQ,tumbleStarts,tumbleLengths,p,opsMovie];*)
(**)
(*(** plot: intensity **)*)
(*opsMovie={AspectRatio->1*)
(*,Epilog->{*)
(*		Text[Style[p["plot","timeLabel"],p["plot","fontSize"]],Scaled[{0.5,-0.18}]]*)
(*		,Rotate[Text[Style[p["plot","intensityLabel"],p["plot","fontSize"]],Scaled[{-0.17,0.5}]],90\[Degree]]*)
(*		}*)
(*,ImagePadding->{ipx,{40,8}}*)
(*,ImageSize->55p["plot","mm"]*)
(*,PlotStyle->AbsoluteThickness[3]*)
(*};*)
(*pI=fig3PlotIntensity[timeData,iData,stepStart,stepEnd,ipx,ops,tTicks,tTicks2,showTumbleShadingQ,tumbleStarts,tumbleLengths,p,opsMovie,1.2];*)
(**)
(*(** plot: intensity spectrum **)*)
(*pSpectrum=siMoviePlotIntensitySpectrum[step,freqData,trajDataAll,ageString,p,ipx,runColor];*)
(**)
(*plot=Grid[{{Grid[{{pRT},{pTheta},{pPhi}}],pI,pSpectrum}}];*)
(**)
(*pthFrames=FileNameJoin[{$HomeDirectory,"Desktop","mmamovie_dynamicPlot"}];*)
(*Export[FileNameJoin[{pthFrames,"p_"<>IntegerString[counter++,10,IntegerLength[nFrames]+1]<>".png"}],plot];*)
(**)
(*,{stepPlot,nFrames}];*)
(*][[1]],"\[ThinSpace]s"];*)


(* ::Input:: *)
(*(** render video with ffmpeg **)*)
(*fps=50;*)
(*renderVideo[pthFrames,FileNameJoin[{pthVideo,"Supplementary Videos","Supplementary Video 9 plots young.mp4"}],fps];*)


(* ::Subsection:: *)
(*Fig. 3 (right column): aged bacterium*)


(* ::Subsubsection::Closed:: *)
(*simulation parameter history*)


(* ::Input:: *)
(*(** 10 um/s - slightly wider; 20 um/s - shorter **)*)
(*{u0Dim,runTime0,tumbleTime0}={5.0*^-6,0.9,0.1};        (** aged bacterium;  w=0.11 **)*)
(*{u0Dim,runTime0,tumbleTime0}={5.0*^-6,0.9,0.3};        (** aged bacterium;  w=0.12 **)*)
(*{u0Dim,runTime0,tumbleTime0}={5.0*^-6,1.2,0.1};        (** aged bacterium;  w=0.09 **)*)
(*{u0Dim,runTime0,tumbleTime0}={5.0*^-6,1.5,0.1};        (** aged bacterium;  w=0.08 **)*)
(*{u0Dim,runTime0,tumbleTime0}={10.0*^-6,1.5,0.1};      (** aged bacterium;  w=0.08 oscillatory intensity **)*)
(*{u0Dim,runTime0,tumbleTime0}={5.0*^-6,1.5,0.05};        (** aged bacterium;  w=0.06 **)*)
(*{u0Dim,runTime0,tumbleTime0}={1.0*^-6,0.2,0.4};         (** aged bacterium; w=0.11 **)*)
(*(*{u0Dim,runTime0,tumbleTime0}={1.0*^-6,0.4,0.4};         (** aged bacterium; w=0.11 **)*)
(*{u0Dim,runTime0,tumbleTime0}={1.0*^-6,0.4,0.2};         (** aged bacterium; w=0.10 **)*)
(*{u0Dim,runTime0,tumbleTime0}={1.0*^-6,0.5,0.2}; *)      (** aged bacterium; w=0.11 **)*)
(*{u0Dim,runTime0,tumbleTime0}={1.0*^-6,0.2504,0.3594};         (** aged bacterium; w=0.11 **)*)
(*{u0Dim,runTime0,tumbleTime0}={1.0*^-6,0.7,0.03};         (** aged bacterium; w=0.11 (great) **)*)
(*{u0Dim,runTime0,tumbleTime0}={1.5*^-6,2.5,0.5};*)
(*{u0Dim,runTime0,tumbleTime0}={4.0*^-6,0.9,0.1};      (** for friction=0.25 friction0, w=0.23 **)*)
(*{u0Dim,runTime0,tumbleTime0}={1.0*^-6,0.9,0.1};      (** for friction=0.25 friction0, w=0.23 **)*)
(*{u0Dim,runTime0,tumbleTime0}={8.0*^-6,0.9,0.1};      (** for friction=0.25 friction0, w=0.20 **)*)
(*{u0Dim,runTime0,tumbleTime0}={5.0*^-6,0.9,0.01};      (** for friction=0.25 friction0, w=0.08 **)*)
(*{u0Dim,runTime0,tumbleTime0}={6.8*^-6,1.0,0.013};      (** for friction=0.195 friction0, w=0.08 **)*)
(*{u0Dim,runTime0,tumbleTime0}={7.2*^-6,0.9,0.011};      (** for friction=0.195 friction0, w=0.08 **)*)
(**)


(* ::Subsubsection:: *)
(*run simulation*)


(* ::Input:: *)
(*{u0Dim,runTime0,tumbleTime0,frictionFactor}={3.2*^-6,1.62,0.01,1.0}(*{3.2*^-6(*7.2*^-6*),0.9,0.011}*);      (** for friction=0.195 friction0, w=0.08 **)*)
(*loadDataQ=False;*)
(*{timeDataOld,freqDataOld,trajDataAllOld}=getSimulationData[False,u0Dim,runTime0,tumbleTime0,frictionFactor,loadDataQ,p];*)


(* ::Subsubsection:: *)
(*save simulation as new default*)


(* ::Input:: *)
(*(** Note: only run this when you are sure you want to keep the results! **)*)
(*saveNewSimulationResultAsDefault["old",p];*)


(* ::Subsubsection:: *)
(*plot data*)


(* ::Input:: *)
(*nIntensityImages=4;*)
(*trajIndex=3;*)
(*{showTumbleShadingQ,colorizeTumblesQ}={True,True};*)
(*\[Theta]PlotMaxDeg=60;*)
(*\[Phi]OffsetHemisphere=0.5\[Pi];*)
(*\[Phi]OffsetTimeseries=0.0\[Pi];*)
(*ageString="old";*)
(*{fwhmPointOld,thetaData,phiData,subplots}=fig3PlotState[ageString,timeDataOld,freqDataOld,trajDataAllOld,showTumbleShadingQ,colorizeTumblesQ,\[Theta]PlotMaxDeg,dt,pthOutFig3,p,nIntensityImages,trajIndex,\[Phi]OffsetHemisphere,\[Phi]OffsetTimeseries];*)
(*fwhmPointOld*)


(* ::Subsubsection:: *)
(*SI Video: old bacterium with plots*)


(* ::Input:: *)
(*(** function from droplet optics **)*)
(*solution=solveHemisphereLC[];*)


(* ::Input:: *)
(*(** create video frames **)*)
(*pthFrames=FileNameJoin[{$HomeDirectory,"Desktop","mmamovie_runTumbleSimulation_"<>ageString}];*)
(*{tmax,tfinalSim}={20,200};*)
(**)
(*Print["movie frame generation took: ",AbsoluteTiming[*)
(*createMovieFramesSimulation[trajectories[[trajIndex,stepsTransient;;stepsTransient+Round[tmax/tfinalSim nSteps]]],ageString,dt,trajDataAllOld,timeDataOld,freqDataOld,pthFrames,p,tmax,2*)
(*];*)
(*][[1]],"\[ThinSpace]s"];*)


(* ::Input:: *)
(*(** render video with ffmpeg **)*)
(*fps=50;*)
(*renderVideo[pthFrames,FileNameJoin[{pthVideo,"Supplementary Videos","Supplementary Video 9 Run Tumble Sphere "<>ageString<>" full.mp4"}],fps];*)


(* ::Subsubsection::Closed:: *)
(*SI Video: old bacterium*)


(* ::Text:: *)
(*$name = 'Supplementary Video 10'*)
(*$fps = '50'*)
(*ffmpeg -y -framerate "$fps" -i "C:\Users\Jan\Dropbox (MIT)\Bacteria Droplet Manuscript\Code\Run Tumble\mmamovie_2\p_%05d.png" -vf scale="trunc(iw/2)*2:trunc(ih/2)*2" -vcodec libx264 -pix_fmt yuv420p -preset slow -crf 10 -r "$fps"  "C:\Users\Jan\Dropbox (MIT)\Bacteria Droplet Manuscript\Videos\Supplementary Videos\$name.mp4"*)


(* ::Input:: *)
(*tmax=20;*)
(*AbsoluteTiming[*)
(*createMovieFramesPosition[trajectories[[trajIndex,stepsTransient;;stepsTransient+Round[tmax/tfinalSim nSteps];;2]],2,dt];*)
(*]*)


(* ::Subsubsection:: *)
(*SI Video: intensity images*)


(* ::Input:: *)
(*(** get intensity images **)*)
(*{timeDuration,tStart,stepStart,stepEnd}=fig3SetTime[dt];*)
(*nIntensityImages=stepEnd-stepStart+1;*)
(*intensityImages0=getIntensityImages[thetaData,phiData,{stepStart,stepEnd},p,nIntensityImages];*)
(**)
(*(** create video frames **)*)
(*pthFrames=FileNameJoin[{$HomeDirectory,"Desktop","mmamovie_intensity_images_old"}];*)
(*If[!DirectoryQ@#,CreateDirectory@#]&@pthFrames;*)
(*counter=0;*)
(*Do[*)
(*	Export[FileNameJoin[{pthFrames,"p_"<>IntegerString[counter++,10,IntegerLength[Length[intensityImages0]]+1]<>".png"}],intensityImages0[[i]]]*)
(*,{i,Length[intensityImages0]}];*)
(**)
(*(** render video with ffmpeg **)*)
(*fps=50;*)
(*renderVideo[pthFrames,FileNameJoin[{pthVideo,"Supplementary Videos","Supplementary Video 9 intensity old.mp4"}],fps];*)


(* ::Chapter:: *)
(*Figure 4*)


(* ::Section:: *)
(*functions and constants*)


(* ::Subsection:: *)
(*fig4*)


(* ::Input:: *)
(*fig4[p];*)


(* ::Subsubsection::Closed:: *)
(*code*)


fig4[p_]:=Module[{freqsEx,ampsEx,widthsEx,meanSpectraEx,meanSpectraSim,meanMap,stdMap,meanMapInterpolation,stdMapInterpolation,plotOpsFig4,pthOutFig4},
	{freqsEx,ampsEx,widthsEx,meanSpectraEx,meanSpectraSim,meanMap,stdMap,meanMapInterpolation,stdMapInterpolation,plotOpsFig4,pthOutFig4}=fig4PrepareData[p];
	
	fig4abPlot[freqsEx,ampsEx,meanSpectraEx,meanSpectraSim,pthOutFig4,plotOpsFig4,p];
	fig4cPlot[widthsEx,plotOpsFig4,pthOutFig4,p];
	fig4dPlot[meanMap,stdMap,meanMapInterpolation,widthsEx,plotOpsFig4,pthOutFig4,p];
	fig4ePlot[meanMapInterpolation,stdMapInterpolation,widthsEx,plotOpsFig4,pthOutFig4,p];
]


(* ::Subsection::Closed:: *)
(*fig4abPlot: Spectra of young and aged bacteria*)


(* ::Input:: *)
(*fig4abPlot[freqsEx,ampsExAmp,meanSpectraEx,meanSpectraSim,pthOut,plotOpsFig4,p]*)


(* ::Subsubsection::Closed:: *)
(*test: plot std data*)


(* ::Text:: *)
(*TODO : integrate into main function*)


(* ::Input:: *)
(*pTicks=Table[{t,DecimalForm[t,{2,1}]},{t,0,1,0.5}];*)
(*pTicks2={#[[1]],""}&/@pTicks;*)
(*yMin=-0.02;*)
(*yMax=4.5;*)
(*fwhmPoints={{0.15,0.5},{0.07,0.5}};*)
(*ops={PlotRange->{{-0.01,1.01}fMax,(*{-0.02,1.05}*)All},Frame->True,AspectRatio->1,FrameLabel->{p["plot","frequencyLabel"],p["plot","spectrumLabel"]},FrameStyle->Directive[Black,18,AbsoluteThickness[2],"Arial"],ImageSize->{Automatic,300},FrameTicks->{{pTicks,pTicks2},Automatic(*{fTicks,fTicks2}*)}};*)
(*dayIndices={1,6};*)
(*{plots,stdplots}=Table[*)
(*	color=If[plotIndex==1,colorYoung,colorOld];*)
(*	meanPSD={freqs,#/Max[#]}\[Transpose]&@meanSpectraEx[[dayIndices[[plotIndex]]]];*)
(*	meanPSDpSTD={freqs,(#[[1]]+#[[2]])/Max[#[[1]]]}\[Transpose]&@{meanSpectraEx[[dayIndices[[plotIndex]]]],stdSpectra[[dayIndices[[plotIndex]]]]};*)
(*	meanPSDmSTD={freqs,(#[[1]]-#[[2]])/Max[#[[1]]]}\[Transpose]&@{meanSpectraEx[[dayIndices[[plotIndex]]]],stdSpectra[[dayIndices[[plotIndex]]]]};*)
(*	meanPSDmSTD=Max[#,yMin]&/@meanPSDmSTD;*)
(*	pMean=ListLinePlot[meanPSD,Evaluate[ops],PlotStyle->color,Epilog->{Red,AbsolutePointSize[10],Point[#],AbsoluteThickness[2],Dashed,Line[{{-0.01,0.5},#}],Line[{#,{#[[1]],yMin}}]}&@fwhmPoints[[plotIndex]]];*)
(*	pMeanStd=ListLinePlot[{meanPSD,meanPSDpSTD,meanPSDmSTD},PlotRange->{{-0.01,1.01}fMax,{yMin,yMax}},Evaluate[ops],PlotStyle->color,Filling->{2->{3}},FillingStyle->Lighter@color,Epilog->{Red,AbsolutePointSize[10],Point[#],AbsoluteThickness[2],Dashed,Line[{{-0.01,0.5},#}],Line[{#,{#[[1]],yMin}}]}&@fwhmPoints[[plotIndex]]];*)
(*{pMean,pMeanStd}*)
(*,{plotIndex,{1,2}}]\[Transpose];*)


(* ::Subsubsection:: *)
(*code*)


fig4abPlot[freqsEx_,ampsExAmp_,meanSpectraEx_,meanSpectraSim_,pthOut_,plotOpsFig4_,p_]:=Module[
{fTicks,fTicks2,sTicks,sTicks2,yMin,yMax,fwhmPointYoung,fwhmPointOld,fwhmPoints,dayIndices,ops4ab,color,meanPSDex,pMeanEx,meanPSDsim,pMeanSim,plots4ab},
	
	fTicks=Table[If[Mod[f,0.5]==0,{f,DecimalForm[f,{3,1}],{0.02,0}},{f,"",{0.01,0}}],{f,0,1,0.1}];
	fTicks2={#[[1]],"",#[[3]]}&/@fTicks;
	sTicks=Table[If[Mod[s,1.0]==0,{s,DecimalForm[s,{2,1}],{0.02,0}},{s,"",{0.01,0}}],{s,0,2.0,0.5}];
	sTicks2={#[[1]],""}&/@sTicks;
	{yMin,yMax}={-0.02,4.5};
	
	{fwhmPointYoung,fwhmPointOld}={{0.23,0.5},{0.08,0.5}};
	fwhmPoints={fwhmPointYoung,fwhmPointOld};
	
	dayIndices={1,6};
	
	plots4ab=Table[
		ops4ab={
			Evaluate[plotOpsFig4]
			,AspectRatio->1
			,FrameTicks->{{sTicks,sTicks2},{fTicks,fTicks2}}
			,PlotRange->{{0,1}p["fourier","fMax"],{-0.05,1.6}}
			,PlotRangeClipping->False  (** needed to set frame labels via Epilog **)
			,Epilog->{
				(** width highlight **)
				{Red,AbsolutePointSize[5],Point[#],AbsoluteThickness[p["plot","lineThickness"]],Dashing[0.02],HalfLine[{#,{0,0.5}}],HalfLine[{#,{#[[1]],0}}]}&@fwhmPoints[[plotIndex]]
				(** frame labels **)
				,{
					Text[Style[p["plot","frequencyLabel"],p["plot","fontSize"]],Scaled[{0.5,-0.3}]]
					,Rotate[Text[Style["Intensity spectrum",p["plot","fontSize"]],Scaled[{-0.37,0.5}]],90\[Degree]]
				}
			}
		};
		color=p["color","youngOld"][[plotIndex]];
		
		(** experiments **)
		meanPSDex={freqsEx,#/(Max[#]ampsExAmp[[dayIndices[[plotIndex]]]])}\[Transpose]&@meanSpectraEx[[dayIndices[[plotIndex]]]];
		pMeanEx=ListLinePlot[meanPSDex,Evaluate[ops4ab],PlotStyle->color];
		
		(** simulations **)
		meanPSDsim=meanSpectraSim[[plotIndex]];
		pMeanSim=ListLinePlot[meanPSDsim,Evaluate[ops4ab],PlotStyle->Directive[Darker[color,0.5],Dashed]];
		
		(** ex+sim comparison plot **)
		Show[{pMeanEx,pMeanSim},ImageSize->45p["plot","mm"]]
	,{plotIndex,{1,2}}];
	
	Print[Grid[{plots4ab}]];
	
	(** export plots **)
	MapThread[Export[FileNameJoin[{pthOut,#1<>".png"}],#2,Background->None,ImageResolution->300]&,{{"figure_4a","figure_4b"},plots4ab}];
	MapThread[Export[FileNameJoin[{pthOut,#1<>".pdf"}],#2,Background->None]&,{{"figure_4a","figure_4b"},plots4ab}];
]


(* ::Subsection::Closed:: *)
(*fig4cPlot: Width as a function of age*)


(* ::Input:: *)
(*fig4cPlot[widthsEx,plotOpsFig4,pthOut,p]*)


(* ::Subsubsection::Closed:: *)
(*check fit*)


(* ::Input:: *)
(*Manipulate[*)
(*	p1=Show[*)
(*		ListLinePlot[{freqsEx,#/Max[#]&@meanSpectraEx[[i]]}\[Transpose],Evaluate[p["plot","ops"]],PlotStyle->Gray,PlotRange->All],*)
(*		Plot[iFits[[i]][f],{f,0,1},PlotRange->All]*)
(*		,PlotRange->{{0,1},All}*)
(*	];*)
(*	p2=Show[*)
(*		ListLinePlot[{freqsEx,#/(Max[#]ampsExAmp[[i]])&@meanSpectraEx[[i]]}\[Transpose],Evaluate[p["plot","ops"]],PlotStyle->Gray,PlotRange->All],*)
(*		Plot[iFitsAmp[[i]][f]/ampsExAmp[[i]],{f,0,1},PlotRange->All,PlotStyle->ColorData[97,2]]*)
(*		,PlotRange->{{0,1},All}*)
(*	];*)
(*	Grid[{#}&/@{p1,p2}]*)
(*,{i,1,p["experiment","nDays"],1,Appearance->"Open"}*)
(*,TrackedSymbols:>{i}]*)


(* ::Subsubsection::Closed:: *)
(*test: plot fit*)


(* ::Input:: *)
(*widthFit=LinearModelFit[ageWidthData,t,t];*)


(* ::Subsubsection:: *)
(*code*)


fig4cPlot[widthsEx_,plotOpsFig4_,pthOut_,p_]:=Module[{aTicks,aTicks2,wTicks,wTicks2,ageWidthData,plot4c},
	
	(** frame ticks **)
	aTicks=Range[0,5];
	aTicks2={#,""}&/@aTicks;
	wTicks=Table[{w,ToString@DecimalForm[w,{2,1}]},{w,0,0.3,0.1}];
	wTicks2={#[[1]],""}&/@wTicks;
	
	ageWidthData={Range[0,5],widthsEx}\[Transpose];
	
	plot4c=ListPlot[ageWidthData
		,Evaluate[plotOpsFig4]
		,AspectRatio->1
		,ColorFunction->(ageColormap[#]&)
		,FrameTicks->{{wTicks,wTicks2},{aTicks,aTicks2}}
		,PlotRange->{{0,5},{0,0.3}}
		,PlotRangeClipping->False  (** needed to set frame labels via Epilog **)
		,PlotRangePadding->{Scaled[0.1],{0,Scaled[0.05]}}
		,PlotStyle->Directive[AbsolutePointSize[4]]
		,Epilog->{
			(** frame labels **)
			Text[Style["Age (days)",p["plot","fontSize"]],Scaled[{0.5,-0.3}]]
			,Rotate[Text[Style[p["plot","widthLabel"],p["plot","fontSize"]],Scaled[{-0.37,0.5}]],90\[Degree]]
		}
	];
	plot4c=Show[plot4c,ImageSize->45 p["plot","mm"],Background->None];  (** for correct image size in export **)
	Print[plot4c];
	
	(** export plot **)
	Export[FileNameJoin[{pthOut,"plot4c.png"}],plot4c,ImageResolution->300];
	Export[FileNameJoin[{pthOut,"plot4c.pdf"}],plot4c];
	
]


(* ::Subsection::Closed:: *)
(*fig4dPlot: Width/Angle Mapping*)


(* ::Input:: *)
(*fig4dPlot[meanMap,stdMap,meanMapInterpolation,widthsEx,plotOpsFig4,pthOutFig4,p];*)


(* ::Subsubsection::Closed:: *)
(*alternative visualization*)


(* ::Input:: *)
(*centerBox[{x0_,y0_},{lx_,ly_}]:=Rectangle[{x0-0.5lx,y0-0.5ly},{x0+0.5lx,y0+0.5ly}]*)
(*ListLinePlot[meanMap*)
(*,AspectRatio->1*)
(*,PlotRange->{{0.06,0.26},{5,26}}*)
(*,PlotStyle->Black*)
(*,Prolog->{RGBColor[0.64, 0.64, 0.64],MapThread[centerBox[#1,#2]&,{meanMap,2stdMap}]}*)
(*,Frame->True*)
(*,FrameLabel->{"Width (Hz)","Angle (\[Degree])"}*)
(*,FrameStyle->Directive[Black,20,AbsoluteThickness[2]]*)
(*]*)


(* ::Subsubsection:: *)
(*code*)


Clear[fig4dPlot];
fig4dPlot[meanMap_,stdMap_,meanMapInterpolation_,widthsEx_,plotOpsFig4_,pthOutFig4_,p_]:=Module[{stdUpMap,stdDownMap,\[Theta]Ticks,\[Theta]Ticks2,wTicks,wTicks2,plot4d},
	
	wTicks=Range[0.05,0.25,0.1];
	wTicks2={#,""}&/@wTicks;
	\[Theta]Ticks=Range[0,30,5];
	\[Theta]Ticks2={#,""}&/@\[Theta]Ticks;
	
	stdUpMap={meanMap[[All,1]],meanMap[[All,2]]+stdMap[[All,2]]}\[Transpose];
	stdDownMap={meanMap[[All,1]],meanMap[[All,2]]-stdMap[[All,2]]}\[Transpose];
	
	plot4d=Show[
		ListLinePlot[Sort/@{meanMap,stdUpMap,stdDownMap}
			,Evaluate[plotOpsFig4]
			,AspectRatio->1
			,Filling->{2->{3}}
			,PlotStyle->Evaluate[Directive[Black,AbsoluteThickness[#]]&/@({1,0.5,0.5}p["plot","lineThickness"])]
			,PlotRange->{{0.045,0.27},{5,26}}
			,PlotRangeClipping->False  (** needed to set frame labels via Epilog **)
			,FrameTicks->{{\[Theta]Ticks,\[Theta]Ticks2},{wTicks,wTicks2}}
			,Epilog->{
				(** frame labels **)
				Text[Style[p["plot","widthLabel"],p["plot","fontSize"]],Scaled[{0.5,-0.3}]]
				,Rotate[Text[Style["Angle \[LeftAngleBracket]\[Theta]\[RightAngleBracket] (\[Degree])",p["plot","fontSize"]],Scaled[{-0.37,0.5}]],90\[Degree]]
			}
		]
		,
		ListPlot[{{#,meanMapInterpolation[#]}}&/@widthsEx,PlotStyle->Evaluate[Directive[#,AbsolutePointSize[6]]&/@p["color","ages"]]]
	];
	plot4d=Show[plot4d,ImageSize->45 p["plot","mm"],Background->None];  (** for correct image size in export **)
	Print[plot4d];
	
	Export[FileNameJoin[{pthOutFig4,"fig4d.png"}],plot4d,ImageResolution->300];
	Export[FileNameJoin[{pthOutFig4,"fig4d.pdf"}],plot4d];
]


(* ::Subsection::Closed:: *)
(*fig4ePlot: Angle as a function of age*)


(* ::Input:: *)
(*fig4ePlot[meanMapInterpolation,stdMapInterpolation,widthsEx,plotOpsFig4,pthOutFig4,p];*)


(* ::Subsubsection::Closed:: *)
(*test: plot all surfaces and intersections*)


(* ::Input:: *)
(*{widthInterpolationLinLinLin,angleInterpolationLinLinLin}=getScanDataInterpolantsLinLinLin[dataScan];*)


(* ::Input:: *)
(*iSave=2;*)
(*linearlinearlinearQ=False;*)
(**)
(*(** {tR, tT, u0} **)*)
(*(** found manually in previous section **)*)
(*parameterPointsStartEnd={{0.1,0.09,10},{2.5,0.5,1.5}};    *)
(*(** found manually in section Fig. 3 **)*)
(*parameterPointsStartEnd={{0.9,0.1,10.5},{0.9,0.01,5}};*)
(*(** found manually to be large enough to encompass all widths **)*)
(*parameterPointsStartEnd={{0.9,0.1,11},{0.9,0.01,4.0}};*)
(*parameterPointsStartEnd={{0.9,0.1,11},{0.9,0.005,1}};*)
(*(** biologically implausible due to increase in absolute tumble time **)*)
(*(*parameterPointsStartEnd={{0.1,0.09,10},{1.1,0.05,1.0}};*)*)
(*(*parameterPointsStartEnd={{1.0,0.1,10},{1.0,0.1,1.0}};*)        (** only speed decreases 10 -> 1 with {1.0,0.1} **)*)
(*(*parameterPointsStartEnd={{0.1,0.09,10},{0.1,0.09,1.0}};*)    (** only speed decreases 10 -> 1 with {0.1,0.09} from intersection **)*)
(**)
(*parameterEvolution=If[linearlinearlinearQ,lineLinearLinearLinear,lineLogLogLinear];*)
(*{widths,meanAngles,parameterPoints}=getWidthAngleParameterPoints[parameterPointsStartEnd,iSave,widthsEx,dataScan,p["color","ages"]];*)
(**)
(*showIntersectionPointsQ=True;*)
(*plot=Show[*)
(*plotWidthSave,*)
(*plotAngleSave,*)
(*Graphics3D[{*)
(*(** end points **)*)
(*{Black,Sphere[{Log10[#[[1]]],Log10[#[[2]]],#[[3]]},{0.05,0.05,1.^3}]&/@(*parameterPoints[[{1,-1}]]*)parameterPointsStartEnd},*)
(*(** intersection points **)*)
(*If[showIntersectionPointsQ,{Black,Sphere[{Log10[#[[1]]],Log10[#[[2]]],#[[3]]},{0.05,0.05,1.^3}]&/@parameterPoints},Nothing],*)
(*(** curve between first and last point **)*)
(*{Black,Line[(*{Log10[#\[LeftDoubleBracket]1\[RightDoubleBracket]],Log10[#\[LeftDoubleBracket]2\[RightDoubleBracket]],#\[LeftDoubleBracket]3\[RightDoubleBracket]}&/@*)(parameterEvolution[#]&/@Subdivide[0.0,1,100])]}*)
(*(** selected point **)*)
(*(*,{Red,Sphere[{Log10[#\[LeftDoubleBracket]1\[RightDoubleBracket]],Log10[#\[LeftDoubleBracket]2\[RightDoubleBracket]],#\[LeftDoubleBracket]3\[RightDoubleBracket]},{0.05,0.05,1.^3}]&@parameterPoints[[iSave]]}*)*)
(*}]*)
(*,AxesLabel->{"\!\(\*SubscriptBox[\(t\), \(Run\)]\) (s)","\!\(\*SubscriptBox[\(t\), \(Tumble\)]\) (s)",Rotate["u (\[Mu]m/s)",\[Pi]/2]}*)
(*,AxesStyle->Directive[Black,11]*)
(*,BoxStyle->Directive[Black]*)
(*,BoundaryStyle->None*)
(*,ImageSize->220*)
(*,Lighting->{White}*)
(*,PlotRange->{Log10@{0.09,9.`},Log10@{0.01,1},{1,uMax}}*)
(*,ScalingFunctions->{"Log10","Log10","Linear"}*)
(*]*)
(**)
(*aTicks=Range[0,5];*)
(*aTicks2={#,""}&/@aTicks;*)
(*wTicks=Table[{w,ToString@DecimalForm[w,{2,1}]},{w,0,0.3,0.1}];*)
(*wTicks2={#[[1]],""}&/@wTicks;*)
(*\[Theta]Ticks=Range[0,60,10];*)
(*\[Theta]Ticks2={#,""}&/@\[Theta]Ticks;*)
(*pDayAngles=ListPlot[meanAngles,Frame->True,FrameLabel->{"age (day)","angle \!\(\**)
(*StyleBox[\"\[LeftAngleBracket]\",\nFontSlant->\"Italic\"]\)\!\(\**)
(*StyleBox[\"\[Theta]\",\nFontSlant->\"Italic\"]\)\!\(\**)
(*StyleBox[\"\[RightAngleBracket]\",\nFontSlant->\"Italic\"]\) (\[Degree])"},PlotRange->{{1,6},{0,All}},PlotRangePadding->{Scaled[.1],{0,Scaled[.23]}},AspectRatio->1,FrameStyle->Directive[Black,AbsoluteThickness[2lineThickness],20],PlotStyle->AbsolutePointSize[10],ImageSize->300,FrameTicks->{{\[Theta]Ticks,\[Theta]Ticks2},{aTicks,aTicks2}},ColorFunction->(ageColormap[#]&)];*)
(*pWidthAngles=ListPlot[MapThread[Style[{#1,#2},#3]&,{widths,meanAngles,p["color","ages"]}],Frame->True,FrameLabel->{"width \!\(\**)
(*StyleBox[SubscriptBox[\"f\", \"HM\"],\nFontSlant->\"Italic\"]\) (Hz)","angle \[LeftAngleBracket]\[Theta]\[RightAngleBracket] (\[Degree])"},PlotRange->{{0.01,All},{0,All}},PlotRangePadding->{{0,Scaled[.15]},{0,Scaled[.23]}},AspectRatio->1,FrameStyle->Directive[Black,AbsoluteThickness[2lineThickness],20],PlotStyle->AbsolutePointSize[10],ImageSize->300,Axes->False,FrameTicks->{{\[Theta]Ticks,\[Theta]Ticks2},{wTicks,wTicks2}}];*)
(*pDayWidths=ListPlot[widths,Frame->True,FrameLabel->{"age (day)","width \!\(\**)
(*StyleBox[SubscriptBox[\"f\", \"HM\"],\nFontSlant->\"Italic\"]\) (Hz)"},PlotRange->{{1,6},{0,All}},PlotRangePadding->{Scaled[.1],{0,Scaled[.23]}},AspectRatio->1,FrameStyle->Directive[Black,AbsoluteThickness[2lineThickness],20],PlotStyle->AbsolutePointSize[10],ImageSize->300,FrameTicks->{{wTicks,wTicks2},{aTicks,aTicks2}},ColorFunction->(ageColormap[#]&)];*)
(**)
(*(*plot=Grid[{{pDayAngles,"   ",pWidthAngles,"   ",pDayWidths}}]*)*)
(*(*plot=Grid[{{pDayWidths,"   ",pDayAngles}}]*)*)
(*plot=Grid[{{pDayWidths,"   ",pWidthAngles,"   ",pDayAngles}}];*)
(**)
(**)
(*plot4x=*)
(*ListPlot[MapThread[Style[{#1,#2},#3]&,{widths,meanAngles,p["color","ages"]}],PlotStyle->Directive[AbsolutePointSize[4]],Frame->True,FrameLabel->{"Width \!\(\**)
(*StyleBox[SubscriptBox[\"f\", \"HM\"],\nFontSlant->\"Italic\"]\) (Hz)","Angle \[LeftAngleBracket]\[Theta]\[RightAngleBracket] (\[Degree])"},PlotRange->{{0.01,0.3},{0,All}},PlotRangePadding->{{0,Scaled[0.05]},{0,Scaled[0.23]}},AspectRatio->1,FrameStyle->Directive[Black,fontSize,AbsoluteThickness[lineThickness]],FrameTicks->{{\[Theta]Ticks,\[Theta]Ticks2},{wTicks,wTicks2}},ImagePadding->ipadSubfiguresFig4,Axes->False,ImageSize->45mm];*)
(**)
(*days=Range[0,5];*)
(*plot4d=*)
(*ListPlot[{days,meanAngles}\[Transpose],PlotStyle->Directive[AbsolutePointSize[4]],ColorFunction->(ageColormap[#]&),*)
(*Frame->True,PlotRange->{{0,5},{0,All}},PlotRangePadding->{Scaled[.1],{0,Scaled[.23]}},AspectRatio->1,FrameStyle->Directive[Black,fontSize,AbsoluteThickness[lineThickness]],FrameLabel->{"Age (days)","Angle \!\(\**)
(*StyleBox[\"\[LeftAngleBracket]\[Theta]\[RightAngleBracket]\",\nFontSlant->\"Italic\"]\) (\[Degree])"},FrameTicks->{{\[Theta]Ticks,\[Theta]Ticks2},{aTicks,aTicks2}},ImagePadding->ipadSubfigures,Axes->False,ImageSize->45mm]*)


(* ::Subsubsection:: *)
(*test: error computation*)


(* ::Text:: *)
(*meanMapInterpolation = \[Theta](w)*)


(* ::Input:: *)
(*Plot[meanMapInterpolation[x],{x,0.0701,0.25}]*)


(* ::Input:: *)
(*(** compute error of angle: \[CapitalDelta]\[Theta] = d\[Theta]/dw \[CapitalDelta]w **)*)
(*angleError=meanMapInterpolation'[widthsEx] stdMapInterpolation[widthsEx]*)


(* ::Subsubsection::Closed:: *)
(*test: visualize start/end point volumes/surfaces*)


(* ::Input:: *)
(*maps=filteredMappings3Inc[[All,3]];*)
(*ListPlot[maps[[;;;;100]]]*)
(*ListPlot[Mean[maps[[;;;;100,All(*,2*)]]]]*)


(* ::Input:: *)
(*ListPlot[filteredMappings[[;;;;100]]]*)


(* ::Input:: *)
(*mappings=filteredMappings;*)
(*Manipulate[*)
(*	ListPlot[mappings[[i]],Frame->True,FrameLabel->{"Width","Mean angle"}]*)
(*,{i,1,Length[mappings],1},TrackedSymbols:>{i}]*)


(* ::Input:: *)
(*(** surfaces **)*)


(* ::Input:: *)
(*ops={ScalingFunctions->{"Log10","Log10",None},Mesh->None};*)
(*plotWidthContours=ContourPlot3D[widthInterpolationLinLinLin[x,y,z],{x,0.09,9},{y,0.01,1},{z,1,10},ContourStyle->Opacity[0.6],Contours->{wMin,wMax},PlotLegends->Table["width = "<>ToString[DecimalForm[w,{2,2}]]<>" Hz",{w,{wMin,wMax}}]]*)


(* ::Input:: *)
(*(** volumes **)*)
(*regionSmallWidth=DiscretizeRegion[ImplicitRegion[widthInterpolationLinLinLin[x, y, z]<wMin\[And](0.09<=x<=9)\[And](0.01<=y<=1)\[And](1<=z<=10),{x,y,z}]];*)
(*regionLargeWidth=DiscretizeRegion[ImplicitRegion[widthInterpolationLinLinLin[x, y, z]>wMax\[And](0.09<=x<=9)\[And](0.01<=y<=1)\[And](1<=z<=10),{x,y,z}]];*)


(* ::Input:: *)
(*RegionPlot3D[{regionSmallWidth,regionLargeWidth}, (*ScalingFunctions -> {"Log10", "Log10", None},*)Mesh->None,PlotStyle->Directive[Opacity[0.5]]]*)


(* ::Input:: *)
(*RegionPlot3D[Evaluate[{#<=Min[widthsEx],#>=Max[widthsEx]}&@widthInterpolationLinLinLin[x, y, z]], {x, 0.09, 9}, {y, 0.01, 1}, {z,*)
(*   1, 10}, ScalingFunctions -> {"Log10", "Log10", None},Mesh->None,PlotStyle->Directive[Opacity[0.5]]]*)


(* ::Subsubsection::Closed:: *)
(*test: point positions relative to intersections of width and angle surfaces*)


(* ::Input:: *)
(*i=1; (** 1 or 2 **)*)
(*goodPlotQualityQ=False;  (** nice but slow; leads to same results **)*)
(**)
(*parameterPointsStartEnd={{0.1,0.09,10},(*{2.5,0.5,1.5},*){1.1,0.05,1.0}};*)
(*meanAngles={ 12.281,3.400};*)
(*{width,meanAngle,parameterPoint}={widthsEx[[{1,-1}]],meanAngles[[{1,-1}]],parameterPointsStartEnd}\[Transpose][[i]]*)
(*widthInterpolation=widthInterpolationLinLinLin@@parameterPoint;*)
(*meanAngleInterpolation=angleInterpolationLinLinLin@@parameterPoint;*)
(*Print["parameterPoint = ",parameterPoint];*)
(*Print["widthInterpolation = ",widthInterpolation];*)
(*Print["meanAngleInterpolation = ",meanAngleInterpolation];*)
(*ops={ScalingFunctions->{"Log10","Log10",None},Mesh->None};*)
(**)
(*(** plots **)*)
(*If[goodPlotQualityQ,*)
(*	plotWidth=ContourPlot3D[widthInterpolationLinLinLin[x,y,z],{x,0.09,9},{y,0.01,1},{z,1,10},Evaluate[ops],ContourStyle->Opacity[0.6],Contours->{width},PlotLegends->{"width = "<>ToString[DecimalForm[width,{2,2}]]<>" Hz"}];*)
(*	plotAngle=ContourPlot3D[angleInterpolationLinLinLin[x,y,z],{x,0.09,9},{y,0.01,1},{z,1,10},Evaluate[ops],ContourStyle->Directive[Opacity[0.6],Red],Contours->{meanAngle},PlotLegends->{"mean angle = "<>ToString[DecimalForm[meanAngle,{3,1}]]<>" \[Degree]"}];*)
(*,*)
(*	plotWidth=ListContourPlot3D[dataScan[[All,{1,2,3,4}]],Evaluate[ops],ContourStyle->Opacity[0.6],Contours->{width},PlotLegends->{"width = "<>ToString[DecimalForm[width,{2,2}]]<>" Hz"}];*)
(*	plotAngle=ListContourPlot3D[dataScan[[All,{1,2,3,5}]],Evaluate[ops],ContourStyle->Directive[Opacity[0.6],Red],Contours->{meanAngle},PlotLegends->{"mean angle = "<>ToString[DecimalForm[meanAngle,{3,1}]]<>" \[Degree]"}];*)
(*];*)
(**)
(*plot=Show[*)
(*	plotWidth*)
(*	,plotAngle*)
(*	,Graphics3D[{*)
(*		(** selected point **)*)
(*		{Darker@Gray,Scale[Sphere[{Log10[#[[1]]],Log10[#[[2]]],#[[3]]},0.05],{1,1,4}](*,{0.05,0.05,0.01}*)&@parameterPoint}*)
(*	}]*)
(*	,PlotRange->{Log10@{0.09,9.`},Log10@{0.01,1},{1,10}}*)
(*	,AxesLabel->{"\!\(\*SubscriptBox[\(t\), \(Run\)]\) (s)","\!\(\*SubscriptBox[\(t\), \(Tumble\)]\) (s)",Rotate["u (\[Mu]m/s)",\[Pi]/2]}*)
(*	,AxesStyle->Directive[Black,12]*)
(*	,BoxStyle->Directive[Black]*)
(*]*)
(*Export[FileNameJoin[{$HomeDirectory,"Desktop","plot_SinglePoint_"<>ToString[i]<>".png"}],plot]*)


(* ::Subsubsection::Closed:: *)
(*test: width(speed)*)


(* ::Input:: *)
(*simData=Import[FileNameJoin[{NotebookDirectory[],"output","output_fig4_scan_test","simAnalysisData_rDroplet_5_3.50_rEcoli_0.5.h5"}],"/Dataset1"];*)


(* ::Input:: *)
(*plotOps={AspectRatio->1,Frame->True,FrameStyle->Directive[Black,fontSize,AbsoluteThickness[lineThickness],"Arial"],*)
(*PlotRangePadding->{Scaled[0.02],{0,Scaled@0.05}},ColorFunction->Function[{t,w},ageColormap[t]],Mesh->All,MeshStyle->AbsolutePointSize[3],ScalingFunctions->{"Reverse",Identity},PlotStyle->AbsoluteThickness[lineThickness],ImagePadding->ipad};*)
(*u0DimData=10^6 GeneralUtilities`AssociationTranspose[simData]["u0Dim"];*)
(*iWidthData=GeneralUtilities`AssociationTranspose[simData]["iWidth"];*)
(**)
(*uTicks=Table[If[Mod[i,2]==0,{i,i,{0.015,0}},{i,"",{0.0075,0}}],{i,1,10,1}];*)
(*uTicks2={#[[1]],"",#[[3]]}&/@uTicks;*)
(*wTicks=Table[If[Mod[w,0.1]==0,{w,DecimalForm[w,{2,1}],{0.015,0}},{w,"",{0.0075,0}}],{w,0,0.3,0.1}];*)
(*wTicks2={#[[1]],"",#[[3]]}&/@wTicks;*)
(**)
(*Show[ListLinePlot[{u0DimData,iWidthData}\[Transpose][[;;10]],Evaluate[plotOps],PlotRange->{{1,10},{0,All}},FrameLabel->{"speed \!\(\**)
(*StyleBox[\"u\",\nFontSlant->\"Italic\"]\) (\!\(\*SuperscriptBox[\(10\), \(-6\)]\) m/s)","width \!\(\**)
(*StyleBox[SubscriptBox[\"f\", \"HM\"],\nFontSlant->\"Italic\"]\) (Hz)"},FrameTicks->{{wTicks,wTicks2},{uTicks,uTicks2}}]*)
(*,ImageSize->45mm]*)


(* ::Subsubsection::Closed:: *)
(*test: width(run time, tumble time)*)


(* ::Input:: *)
(*runtimeData=GeneralUtilities`AssociationTranspose[simData]["tRun"];*)
(*tumbletimeData=GeneralUtilities`AssociationTranspose[simData]["tTumble"];*)
(*iWidthData=GeneralUtilities`AssociationTranspose[simData]["iWidth"];*)
(*ListDensityPlot[{runtimeData,tumbletimeData,iWidthData}\[Transpose],PlotRange->All,InterpolationOrder->0,ColorFunction->GrayLevel,FrameLabel->{"runtime \!\(\*SubscriptBox[\(t\), \(R\)]\) (s)","tumbletime \!\(\*SubscriptBox[\(t\), \(t\)]\) (s)"},PlotLegends->BarLegend[Automatic,LegendLabel->"width (Hz)"],PlotRangePadding->None,FrameStyle->Directive[Black,20,AbsoluteThickness[2]],ScalingFunctions->{"Log","Log"}]*)


(* ::Subsubsection:: *)
(*code*)


Clear[fig4ePlot];
fig4ePlot[meanMapInterpolation_,stdMapInterpolation_,widthsEx_,plotOpsFig4_,pthOutFig4_,p_]:=Module[{angleError,ageWidthData,aTicks,aTicks2,\[Theta]Ticks,\[Theta]Ticks2,plot4e},
	
	(** "error" of angle = computed standard deviation of model curves **)
	angleError=stdMapInterpolation[widthsEx];
	ageWidthData={#}&/@({Range[0,5],MapThread[Around[#1,#2]&,{meanMapInterpolation[widthsEx],angleError}]}\[Transpose]);
	
	(** frame ticks **)
	aTicks=Range[0,5];
	aTicks2={#,""}&/@aTicks;
	\[Theta]Ticks=Range[0,30,5];
	\[Theta]Ticks2={#,""}&/@\[Theta]Ticks;
	
	plot4e=ListPlot[ageWidthData
		,Evaluate[plotOpsFig4]
		,AspectRatio->2
		,PlotRange->{0,All}
		,PlotRangeClipping->False  (** needed to set frame labels via Epilog **)
		,PlotRangePadding->{Scaled[0.08],{Scaled[0.0],Scaled[0.05]}}
		,PlotStyle->Evaluate[Directive[#,AbsolutePointSize[6]]&/@p["color","ages"]]
		,FrameTicks->{{\[Theta]Ticks,\[Theta]Ticks2},{aTicks,aTicks2}}
		,IntervalMarkersStyle->AbsoluteThickness[p["plot","lineThickness"]]
		,Epilog->{
			(** frame labels **)
			Text[Style["Age (days)",p["plot","fontSize"]],Scaled[{0.5,-0.15}]]
			,Rotate[Text[Style["Angle \[LeftAngleBracket]\[Theta]\[RightAngleBracket] (\[Degree])",p["plot","fontSize"]],Scaled[{-0.37,0.5}]],90\[Degree]]
		}
	];
	plot4e=Show[plot4e,ImageSize->45 p["plot","mm"],Background->None];  (** for correct image size in export **)
	Print[plot4e];
	
	Export[FileNameJoin[{pthOutFig4,"fig4e.png"}],plot4e,ImageResolution->300];
	Export[FileNameJoin[{pthOutFig4,"fig4e.pdf"}],plot4e];
]


(* ::Subsection:: *)
(*fig4PrepareData*)


(* ::Input:: *)
(*{freqsEx,ampsEx,widthsEx,meanSpectraEx,meanSpectraSim,meanMap,stdMap,meanMapInterpolation,stdMapInterpolation,plotOpsFig4,pthOutFig4}=fig4PrepareData[p];*)


(* ::Text:: *)
(*Output spectra are cut to fMax = 1.0 to have less issues with plotting later*)


(* ::Subsubsection::Closed:: *)
(*get parameters for Fig. 3 columns*)


(* ::Input:: *)
(*widthInterpolationLogLogLin[Log10@0.9,Log10@0.1, 7.2]*)
(*widthInterpolationLogLogLin[Log10@0.9,Log10@0.011, 7.2]*)


(* ::Subsubsection::Closed:: *)
(*test: Lorentz fits with different free parameters to get width*)


(* ::Text:: *)
(*Theory says that spectrum should not have an offset/shift. Max of noisy spectrum is not guaranteed to be max of true spectrum. Free parameters must be amplitude and width.*)


(* ::Input:: *)
(*{allIntensities,freqsEx,meanSpectraEx,stdSpectraEx,allSpectraEx,weightedMeanSpectraEx,allWeightedSpectraEx}=calculateDailyMeanSpectra[p];*)


(* ::Input:: *)
(**)


(* ::Input:: *)
(*Clear[\[Omega],\[Omega]0,aa,bb,cc,i0];*)
(*fMaxIndex=FirstPosition[freqsEx,_?(#>=p["fourier","fMax"]&)][[1]];*)
(*{indexStart,indexEnd}={2,fMaxIndex};*)
(**)
(*(** fitting normalized average spectrum in intensity **)*)
(*(** free parameters: width, amplitude, shift **)*)
(*iFitsAmpShift=NonlinearModelFit[{freqsEx,#/Max[#]}\[Transpose][[indexStart;;indexEnd]],{aa/(1+(\[Omega]-\[Omega]0)^2/bb^2),aa>0,bb>0,\[Omega]0>=0},{{aa,1},{bb,0.1},{\[Omega]0,0.001}},\[Omega],MaxIterations->1000]&/@weightedMeanSpectraEx;*)
(*Print["widths (free: width,shift,amplitude): ",#["BestFitParameters"][[2,2]]&/@iFitsAmpShift];*)
(**)
(*(** free parameters: width, amplitude **)*)
(*iFitsAmp=NonlinearModelFit[{freqsEx,#/Max[#]}\[Transpose][[indexStart;;indexEnd]],{aa/(1+(\[Omega])^2/bb^2),aa>0,bb>0},{{aa,1},{bb,0.1}},\[Omega],MaxIterations->1000]&/@weightedMeanSpectraEx;*)
(*{ampsExAmp,widthsExAmp}=Transpose[#["BestFitParameters"][[All,2]]&/@iFitsAmp];*)
(*Print["widths (free: width,amplitude): ",#["BestFitParameters"][[2,2]]&/@iFitsAmp];*)
(**)
(*(** free parameters: width, shift **)*)
(*iFitsShift=NonlinearModelFit[{freqsEx,#/Max[#]}\[Transpose][[indexStart;;indexEnd]],{1/(1+(\[Omega]-\[Omega]0)^2/bb^2),bb>0,\[Omega]0>=0},{{bb,0.1},{\[Omega]0,0.001}},\[Omega],MaxIterations->1000]&/@weightedMeanSpectraEx;*)
(*Print["widths (free: width,shift): ",#["BestFitParameters"][[1,2]]&/@iFitsShift];*)
(**)
(*(** free parameters: width **)*)
(*iFits=NonlinearModelFit[{freqsEx,#/Max[#]}\[Transpose][[indexStart;;indexEnd]],{1/(1+(\[Omega])^2/bb^2),bb>0},{{bb,0.1}},\[Omega],MaxIterations->1000]&/@weightedMeanSpectraEx;*)
(*Print["widths (free: width): ",#["BestFitParameters"][[1,2]]&/@iFits];*)


(* ::Input:: *)
(*Show[*)
(*	ListLinePlot[Evaluate[{freqsEx,#/Max[#]}\[Transpose]&/@weightedMeanSpectraEx],PlotRange->{{0,1},All}],*)
(*	Plot[Evaluate[#[f]&/@iFitsAmp],{f,0,1}]*)
(*]*)


(* ::Input:: *)
(*ListLinePlot[{Range[p["experiment","nDays"]],widthsExAmp}\[Transpose],Frame->True,AspectRatio->1,ImageSize->200,PlotRange->{{1,6},{0,All}},PlotRangePadding->{Scaled[.05],{0,Scaled[.05]}},Mesh->All]*)


(* ::Input:: *)
(*#["ParameterTable"]&/@iFitsAmp*)
(*#["ParameterErrors"]&/@iFitsAmp*)


(* ::Subsubsection:: *)
(*test: visualize start/end point volumes/surfaces*)


(* ::Text:: *)
(*linear-linear-linear interpolation will result in the same contours, with the same logarithmized x,y coordinates*)


(* ::Input:: *)
(*dataScan=loadSimulationScanData[p];*)


(* ::Input:: *)
(*zmax=6.5(*10*);*)
(*{widthInterpolationLogLogLin,angleInterpolationLogLogLin}=getScanDataInterpolantsLogLogLin[dataScan];*)
(*{wMin,wMax}={Floor[#[[1]],0.01],Ceiling[#[[2]],0.01]}&@MinMax[widthsEx];*)
(*plotWidthContourMinWidth=ContourPlot3D[widthInterpolationLogLogLin[x,y,z],{x,Log10@0.01,Log10@10},{y,Log10@0.01,Log10@1},{z,1,zmax},Contours->{wMin}];*)
(*plotWidthContourMaxWidth=ContourPlot3D[widthInterpolationLogLogLin[x,y,z],{x,Log10@0.01,Log10@10},{y,Log10@0.01,Log10@1},{z,1,zmax},Contours->{wMax}];*)
(**)
(*dn=100;*)
(*meshPointsMinWidth=First[Cases[plotWidthContourMinWidth, GraphicsComplex[points_, ___] :> points, Infinity]][[;;;;dn]];*)
(*meshPointsMaxWidth=First[Cases[plotWidthContourMaxWidth, GraphicsComplex[points_, ___] :> points, Infinity]][[;;;;dn]];*)
(*nCombos=Length[meshPointsMinWidth]Length[meshPointsMaxWidth]*)
(**)
(*{trTicksLog,ttTicksLog}=setRunTumbleTicks[True];*)
(*plotWidthPointsAndContours=Show[*)
(*ListPointPlot3D[{meshPointsMaxWidth,meshPointsMinWidth},PlotRange->{Log10@{0.0099,10.1},Log10@{0.0099,0.1 1.01},{0.99,zmax 1.001}},Ticks->{trTicksLog,ttTicksLog,{1,5,10}},BoxRatios->{1,1,1}]*)
(*,ContourPlot3D[widthInterpolationLogLogLin[x,y,z],{x,Log10@0.01,Log10@10},{y,Log10@0.01,Log10@0.1},{z,1,zmax},ContourStyle->Opacity[0.6],Contours->{wMin,wMax},PlotLegends->Table["Width = "<>ToString[DecimalForm[w,{2,2}]]<>" Hz",{w,{wMin,wMax}}],Mesh->None]*)
(*]*)


(* ::Subsubsection::Closed:: *)
(*test: get all mappings*)


(* ::Input:: *)
(*Print[AbsoluteTiming[*)
(*		allMappings=Flatten[Table[getWidthAngleMapping[{pStart,pEnd}],{pStart,meshPointsMaxWidth},{pEnd,meshPointsMinWidth}],1]*)
(*	][[1]]];*)


allMappings//Length


(* ::Subsubsection:: *)
(*visualize increasing and decreasing maps*)


(* ::Input:: *)
(*filteredMappingsDec=Select[filteredMappings,monotonousDecreaseQ[#[[3,All,2]]]&];*)
(*filteredMappingsInc=Select[filteredMappings,monotonousIncreaseQ[#[[3,All,2]]]&];*)
(*Print["# Decreasing/Increasing Mappings: ", Length/@{filteredMappingsDec,filteredMappingsInc}];*)


(* ::Input:: *)
(*Graphics3D[{Opacity[0.2]*)
(*,Green,Line[#[[{1,2}]]]&/@filteredMappingsDec*)
(*,Red,Line[#[[{1,2}]]]&/@filteredMappingsInc*)
(*}*)
(*,Axes->True*)
(*,PlotRange->{Log10@{0.0099,10.1},Log10@{0.0099,0.1 1.01},{0.99,zmax 1.001}},(*Ticks->{trTicksLog,ttTicksLog,{1,5,10}},*)BoxRatios->{1,1,1}*)
(*]*)


(* ::Subsubsection::Closed:: *)
(*debug*)


(* ::Input:: *)
(*(** get experimental spectra for all days **)*)
(*{allIntensities,freqsEx,meanSpectraEx,stdSpectraEx,allSpectraEx,weightedMeanSpectraEx,weightedStdSpectraEx,allWeightedSpectraEx}=calculateDailyMeanSpectra[p];*)
(*	*)
(*(** save experimental spectra for all days **)*)
(*saveDataQ=False;*)
(*If[saveDataQ,Do[Export[FileNameJoin[{p["path","intensitySpectraExperiment"],"exSpectraWeighted_"<>ToString[#]<>".h5"}],{freqsEx,allWeightedSpectraEx[[#]]},{"Datasets",{"freqsEx","allWeightedSpectra"}}]&@i,{i,6}];*)
(*];*)
(**)
(*(** fit all noisy normalized experimental spectra with a Lorentzian of variable width and amplitude **)*)
(*Clear[\[Omega],\[Omega]0,aa,bb,cc,i0];*)
(*fMaxIndex=FirstPosition[freqsEx,_?(#>=p["fourier","fMax"]&)][[1]];*)
(*{indexStart,indexEnd}={2,fMaxIndex};*)
(*iFitsAmp=NonlinearModelFit[{freqsEx,#/Max[#]}\[Transpose][[indexStart;;indexEnd]],{aa/(1+(\[Omega])^2/bb^2),aa>0,bb>0},{{aa,1},{bb,0.1}},\[Omega],MaxIterations->1000]&/@weightedMeanSpectraEx;*)
(*{ampsExAmp,widthsExAmp}=Transpose[#["BestFitParameters"][[All,2]]&/@iFitsAmp];*)
(*If[saveDataQ,*)
(*	Export[FileNameJoin[{p["path","widthsExperiment"],"experimental_widths.dat"}],widthsExAmp];*)
(*];*)
(**)
(*(** get numerical spectra for days 1 (young) and 6 (old) **)*)
(*meanSpectraSim=Import[FileNameJoin[{p["path","intensitySpectraNumerical"],#}],"/Dataset1"]&/@{"meanIntensitySpectrumYoung.hdf5","meanIntensitySpectrumOld.hdf5"};*)
(**)
(*(** get sampled contours **)*)
(*{meshPointsMinWidth,meshPointsMaxWidth,{wMin,wMax}}=getSampledWidthAngleContours[widthsExAmp,p];*)


(* ::Input:: *)
(*{meanMap,stdMap,meanMapInterpolation,stdMapInterpolation}=getWidthAngleMappings[Evaluate@meshPointsMinWidth,Evaluate@meshPointsMaxWidth, wMax];*)


(* ::Subsubsection:: *)
(*code*)


fig4PrepareData[p_]:=Module[{allIntensities,freqsEx,meanSpectraEx,stdSpectraEx,allSpectraEx,weightedMeanSpectraEx,weightedStdSpectraEx,allWeightedSpectraEx,
plotOpsFig4,fMaxIndex,indexStart,indexEnd,iFitsAmp,ampsExAmp,widthsExAmp,meanSpectraSim,dataScan,saveDataQ,meshPointsMinWidth,meshPointsMaxWidth,wMin,wMax,
meanMap,stdMap,meanMapInterpolation,stdMapInterpolation,pthOutFig4,widthInterpolationLogLogLin,angleInterpolationLogLogLin},
	
	saveDataQ=False;
	
	(** get experimental spectra for all days **)
	{allIntensities,freqsEx,meanSpectraEx,stdSpectraEx,allSpectraEx,weightedMeanSpectraEx,weightedStdSpectraEx,allWeightedSpectraEx}=calculateDailyMeanSpectra[p];
	
	(** save experimental spectra for all days **)
	If[saveDataQ,
		Do[Export[FileNameJoin[{p["path","intensitySpectraExperiment"],"exSpectraWeighted_"<>ToString[#]<>".h5"}],{freqsEx,allWeightedSpectraEx[[#]]},{"Datasets",{"freqsEx","allWeightedSpectra"}}]&@i,{i,6}];
	];
	
	(** fit all noisy normalized experimental spectra with a Lorentzian of variable width and amplitude **)
	Clear[\[Omega],\[Omega]0,aa,bb,cc,i0];
	fMaxIndex=FirstPosition[freqsEx,_?(#>=p["fourier","fMax"]&)][[1]];
	{indexStart,indexEnd}={2,fMaxIndex};
	iFitsAmp=NonlinearModelFit[{freqsEx,#/Max[#]}\[Transpose][[indexStart;;indexEnd]],{aa/(1+(\[Omega])^2/bb^2),aa>0,bb>0},{{aa,1},{bb,0.1}},\[Omega],MaxIterations->1000]&/@weightedMeanSpectraEx;
	{ampsExAmp,widthsExAmp}=Transpose[#["BestFitParameters"][[All,2]]&/@iFitsAmp];
	
	(** save experimental fitted widths for all days **)
	If[saveDataQ,
		Export[FileNameJoin[{p["path","widthsExperiment"],"experimental_widths.dat"}],widthsExAmp];
	];
	
	(** get numerical spectra for days 1 (young) and 6 (old) **)
	meanSpectraSim=Import[FileNameJoin[{p["path","intensitySpectraNumerical"],#}],"/Dataset1"]&/@{"meanIntensitySpectrumYoung.hdf5","meanIntensitySpectrumOld.hdf5"};
	
	(** get sampled contours **)
	(**TODO: give interpolation its own function **)
	{dataScan,meshPointsMinWidth,meshPointsMaxWidth,{wMin,wMax},widthInterpolationLogLogLin,angleInterpolationLogLogLin}=getSampledWidthAngleContours[widthsExAmp,p];
	{meanMap,stdMap,meanMapInterpolation,stdMapInterpolation}=getWidthAngleMappings[meshPointsMinWidth, meshPointsMaxWidth,wMax,widthInterpolationLogLogLin,angleInterpolationLogLogLin][[2;;]];
	
	(** plot options **)
	plotOpsFig4={
		Evaluate[p["plot","ops"]]
		,ImagePadding->{{60,11},{45,10}}
	};
	
	(** output directory **)
	pthOutFig4=FileNameJoin[{p["path","figs"],"figure 4"}];
	If[!DirectoryQ@#,CreateDirectory@#]&@pthOutFig4;
	
	(** output **)
	{freqsEx[[;;indexEnd]],ampsExAmp,widthsExAmp,weightedMeanSpectraEx[[All,;;indexEnd]],meanSpectraSim,meanMap,stdMap,meanMapInterpolation,stdMapInterpolation,plotOpsFig4,pthOutFig4}
]


(* ::Subsection:: *)
(*getSampledWidthAngleContours*)


(* ::Input:: *)
(*{dataScan,meshPointsMinWidth,meshPointsMaxWidth,{wMin,wMax},widthInterpolationLogLogLin,angleInterpolationLogLogLin}=getSampledWidthAngleContours[widthsEx,p];*)


(* ::Subsubsection::Closed:: *)
(*debug*)


(* ::Input:: *)
(*dataScan=loadSimulationScanData[p];*)


(* ::Input:: *)
(*widthsEx=widthsExAmp;*)
(**)
(*(** get discretized contours **)*)
(*{widthInterpolationLogLogLin2,angleInterpolationLogLogLin2}=getScanDataInterpolantsLogLogLin[dataScan];*)
(*{wMin2,wMax2}={Floor[#[[1]],0.01],Ceiling[#[[2]],0.01]}&@MinMax[widthsEx]*)
(*plotWidthContourMinWidth2=ContourPlot3D[widthInterpolationLogLogLin2[x,y,z],{x,Log10@0.1,Log10@10},{y,Log10@0.01,Log10@1},{z,1,10},Contours->{wMin2}];*)
(*plotWidthContourMaxWidth2=ContourPlot3D[widthInterpolationLogLogLin2[x,y,z],{x,Log10@0.1,Log10@10},{y,Log10@0.01,Log10@1},{z,1,10},Contours->{wMax2}];*)
(**)
(*(** TODO: better sampling of contours than adhoc **)*)
(*dn=100;*)
(*meshPointsMinWidth2=First[Cases[plotWidthContourMinWidth2, GraphicsComplex[points_, ___] :> points, Infinity]][[;;;;dn]];*)
(*meshPointsMaxWidth2=First[Cases[plotWidthContourMaxWidth2, GraphicsComplex[points_, ___] :> points, Infinity]][[;;;;dn]];*)


(* ::Subsubsection:: *)
(*code*)


getSampledWidthAngleContours[widthsEx_,p_]:=Module[
{dataScan,widthInterpolationLogLogLin,angleInterpolationLogLogLin,wMin,wMax,plotWidthContourMinWidth,plotWidthContourMaxWidth,dn,meshPointsMinWidth,meshPointsMaxWidth},
	
	(** get simulation scan data **)
	dataScan=loadSimulationScanData[p];
	
	(** get discretized contours **)
	{widthInterpolationLogLogLin,angleInterpolationLogLogLin}=getScanDataInterpolantsLogLogLin[dataScan];
	{wMin,wMax}={Floor[#[[1]],0.01],Ceiling[#[[2]],0.01]}&@MinMax[widthsEx];
	plotWidthContourMinWidth=ContourPlot3D[widthInterpolationLogLogLin[x,y,z],{x,Log10@0.1,Log10@10},{y,Log10@0.01,Log10@1},{z,1,10},Contours->{wMin}];
	plotWidthContourMaxWidth=ContourPlot3D[widthInterpolationLogLogLin[x,y,z],{x,Log10@0.1,Log10@10},{y,Log10@0.01,Log10@1},{z,1,10},Contours->{wMax}];
	
	(** TODO: better sampling of contours than adhoc **)
	dn=100;
	meshPointsMinWidth=First[Cases[plotWidthContourMinWidth, GraphicsComplex[points_, ___] :> points, Infinity]][[;;;;dn]];
	meshPointsMaxWidth=First[Cases[plotWidthContourMaxWidth, GraphicsComplex[points_, ___] :> points, Infinity]][[;;;;dn]];
	
	(** output **)
	{dataScan,meshPointsMinWidth,meshPointsMaxWidth,{wMin,wMax},widthInterpolationLogLogLin,angleInterpolationLogLogLin}
]


(* ::Subsection::Closed:: *)
(*getScanDataInterpolantsLinLinLin*)


(* ::Input:: *)
(*{widthInterpolationLinLinLin,angleInterpolationLinLinLin}=getScanDataInterpolantsLinLinLin[dataScan];*)


getScanDataInterpolantsLinLinLin[dataScan_]:=Table[Interpolation[dataScan[[All,{1,2,3,index}]],InterpolationOrder->1],{index,{4,5}}]


(* ::Subsection::Closed:: *)
(*getScanDataInterpolantsLogLogLin*)


(* ::Input:: *)
(*{widthInterpolationLogLogLin,angleInterpolationLogLogLin}=getScanDataInterpolantsLogLogLin[dataScan];*)


getScanDataInterpolantsLogLogLin[dataScan_]:=Table[Interpolation[{Log10@#[[1]],Log10@#[[2]],#[[3]],#[[4]]}&/@dataScan[[All,{1,2,3,index}]],InterpolationOrder->1],{index,{4,5}}]


(* ::Subsection::Closed:: *)
(*getWidthAngleMappings*)


(* ::Input:: *)
(*{filteredMappingsDec,meanMap,stdMap,meanMapInterpolation,stdMapInterpolation}=getWidthAngleMappings[meshPointsMinWidth, meshPointsMaxWidth, wMax,widthInterpolationLogLogLin,angleInterpolationLogLogLin];*)


(* ::Subsubsection::Closed:: *)
(*debug*)


(* ::Code:: *)
(*tTumbleMax=0.1;*)
(*Print[AbsoluteTiming[*)
(*	filteredMappings2=Reap[Do[*)
(*	If[(** conditions: initial+final runtimes are larger than tumbletimes + tumbletime decreases + speed decreases **)*)
(*		pStart[[2]]<=Log10[tTumbleMax]\[And]pStart[[1]]>pStart[[2]]\[And]pEnd[[1]]>pEnd[[2]]\[And]pEnd[[2]]<pStart[[2]]\[And]pEnd[[3]]<pStart[[3]]\[And]pStart[[3]]<7.0, *)
(*		mapping=getWidthAngleMapping[{pStart,pEnd},widthInterpolationLogLogLin,angleInterpolationLogLogLin];*)
(*		If[AllTrue[mapping[[All,1]],#<1.01wMax&]\[And]monotonousDecreaseQ[mapping[[All,1]]],*)
(*			Sow[{pStart,pEnd,mapping}]*)
(*		]*)
(*	]*)
(*	,{pStart,meshPointsMaxWidth},{pEnd,meshPointsMinWidth}]][[2,1]];*)
(*][[1]],"\[ThinSpace]s"];*)


widthsEx


Length[filteredMappings2]
ListPlot[filteredMappings2[[;;;;10,3]],PlotRange->All,Frame->True]


Length[filteredMappings]
ListPlot[filteredMappings[[;;10,3]],PlotRange->All,Frame->True]


	Print["# all mappings: ",Length[filteredMappings]];
	filteredMappingsDec=Select[filteredMappings2,monotonousDecreaseQ[#[[3,All,2]]]&];
	filteredMappingsInc=Select[filteredMappings2,monotonousIncreaseQ[#[[3,All,2]]]&];
	Print["# Decreasing/Increasing Mappings: ", Length/@{filteredMappingsDec,filteredMappingsInc}];


filteredMappingsInc


ListPlot[filteredMappingsInc[[;;,3]],PlotRange->All,Frame->True]


{wMin,wMax}


(* ::Input:: *)
(*widthsEx=Flatten[Import[FileNameJoin[{p["path","widthsExperiment"],"experimental_widths.dat"}]]];*)
(*{dataScan,meshPointsMinWidth,meshPointsMaxWidth,{wMin,wMax},widthInterpolationLogLogLin,angleInterpolationLogLogLin}=getSampledWidthAngleContours[widthsEx,p];*)
(*{filteredMappingsDec,meanMap,stdMap,meanMapInterpolation,stdMapInterpolation}=getWidthAngleMappings[meshPointsMinWidth, meshPointsMaxWidth,wMax,widthInterpolationLogLogLin,angleInterpolationLogLogLin];*)
(*filteredMappingsDec//Dimensions*)


(* ::Input:: *)
(*Graphics3D[{Opacity[0.2]*)
(*,Green,Line[#[[{1,2}]]]&/@filteredMappingsDec*)
(*(*,Red,Line[#[[{1,2}]]]&/@filteredMapping3Inc*)*)
(*}*)
(*,Axes->True*)
(*,PlotRange->{Log10@{0.099,10.1},Log10@{0.0099,0.1 1.01},{0.99,10.01}},(*Ticks->{trTicksLog,ttTicksLog,{1,5,10}},*)BoxRatios->{1,1,1}*)
(*]*)


(* ::Input:: *)
(*ListLinePlot[filteredMappingsDec[[All,3]][[;;;;100]],Evaluate[plotOpsFig4],PlotStyle->Gray,AspectRatio->1]*)


(* ::Subsubsection:: *)
(*code*)


getWidthAngleMappings[meshPointsMinWidth_,meshPointsMaxWidth_,wMax_,widthInterpolationLogLogLin_,angleInterpolationLogLogLin_]:=Module[
{tTumbleMax,filteredMappings,mapping,filteredMappingsDec,filteredMappingsInc,maps,meanMap,stdMap,meanMapInterpolation,stdMapInterpolation},
	
	(** get filtered mappings **)
	tTumbleMax=0.1;
	Print[AbsoluteTiming[
		filteredMappings=Reap[Do[
		If[(** conditions: initial+final runtimes are larger than tumbletimes + tumbletime decreases + speed decreases **)
			pStart[[2]]<=Log10[tTumbleMax]\[And]pStart[[1]]>pStart[[2]]\[And]pEnd[[1]]>pEnd[[2]]\[And]pEnd[[2]]<pStart[[2]]\[And]pEnd[[3]]<pStart[[3]]\[And]pStart[[3]]<7.0, 
			mapping=getWidthAngleMapping[{pStart,pEnd},widthInterpolationLogLogLin,angleInterpolationLogLogLin];
			If[AllTrue[mapping[[All,1]],#<1.01wMax&]\[And]monotonousDecreaseQ[mapping[[All,1]]],
				Sow[{pStart,pEnd,mapping}]
			]
		]
		,{pStart,meshPointsMaxWidth},{pEnd,meshPointsMinWidth}]][[2,1]];
	][[1]],"\[ThinSpace]s"];	
	
	Print["# all mappings: ",Length[filteredMappings]];
	filteredMappingsDec=Select[filteredMappings,monotonousDecreaseQ[#[[3,All,2]]]&];
	filteredMappingsInc=Select[filteredMappings,monotonousIncreaseQ[#[[3,All,2]]]&];
	Print["# Decreasing/Increasing Mappings: ", Length/@{filteredMappingsDec,filteredMappingsInc}];
	Print["Error: No decreasing mappings were found!"];
	
	maps=filteredMappingsDec[[All,3]];
	meanMap=Mean[maps];
	stdMap=StandardDeviation[maps];
	meanMapInterpolation=Interpolation[meanMap];
	stdMapInterpolation=Interpolation[{meanMap[[All,1]],stdMap[[All,2]]}\[Transpose]];
	
	(** output **)
	{filteredMappingsDec,meanMap,stdMap,meanMapInterpolation,stdMapInterpolation}
]


(* ::Subsection::Closed:: *)
(*getWidthAngleParameterPoints*)


(* ::Input:: *)
(*{widths,meanAngles,parameterPoints}=getWidthAngleParameterPoints[parameterPointsStartEnd,iSave,widthsEx,cubeData,p];*)


(* ::Subsubsection:: *)
(*debug*)


(* ::Input:: *)
(*i=1;*)
(*parameterPoint=If[i==1,parameterPointsStartEnd[[1]],parameterPointsStartEnd[[-1]]]*)
(*width=widthInterpolationLogLogLin@@{Log10[#[[1]]],Log10[#[[2]]],#[[3]]}&@parameterPoint*)


(* ::Subsubsection::Closed:: *)
(*code*)


getWidthAngleParameterPoints[parameterPointsStartEnd_,widthsEx_,cubeData_,p_]:=Module[
{widthEx,plotWidth,contourMeshRegion,distanceFunction,intersectionPointLogLogLinear,parameterPoint,width,
meanAngleInterpolation,meanAngle,plotAngle,ops,distance,tIntersection},

(** output = loop over all widths **)
(** find intersection points between isowidth surface and evolution curve **)
Table[
	(** get width isosurface and save plot for later **)
	widthEx=widthsEx[[i]];

	(** plot isosurface of Subscript[width, i] **)
	plotWidth=ListContourPlot3D[cubeData[[All,{1,2,3,4}]]
		,BoundaryStyle->None
		,Contours->{#}
		,ContourStyle->Directive[p["color","ages"][[i]],Opacity[0.3]]
		,Mesh->None
		,PlotLegends->{"width = "<>ToString[DecimalForm[#,{5,3}]]<>" Hz"}
		,ScalingFunctions->{"Log10","Log10","Linear"}
	]&@widthEx;

	(** check if experimental width is in simulated data **)
	If[(Min[#]<widthEx<Max[#])&@cubeData[[All,4]],
		(** get width iso-surface **)
		contourMeshRegion=DiscretizeGraphics[Normal[plotWidth/. (Lighting->_):>Lighting->Automatic],Axes->True,Boxed->True];
		
		(** get intersection point **)
		distanceFunction=RegionDistance[contourMeshRegion];
		{distance,tIntersection}=NMinimize[{distanceFunction[parameterEvolution[t]],0<t<1},t];
		If[distance>0.01,Print["Warning (i="<>ToString[i]<>"): Intersection does not exist! Computing closest point on curve to surface"];];
		
		(** assume intersection exists **)
		intersectionPointLogLogLinear=parameterEvolution[t]/.tIntersection[[1]];
		parameterPoint={10^#[[1]],10^#[[2]],#[[3]]}&@intersectionPointLogLogLinear;
		width=widthInterpolationLogLogLin@@{Log10[#[[1]]],Log10[#[[2]]],#[[3]]}&@parameterPoint;
		meanAngleInterpolation=angleInterpolationLogLogLin@@{Log10[#[[1]]],Log10[#[[2]]],#[[3]]}&@parameterPoint;
		meanAngle=meanAngleInterpolation;
		
		plotAngle=ListContourPlot3D[cubeData[[All,{1,2,3,5}]]
			,BoundaryStyle->None
			,Contours->{#}
			,ContourStyle->Directive[Opacity[0.3],ColorData["SunsetColors"][(0.9-0.0)/(55-5)meanAngleInterpolation]]
			,Lighting->{White}
			,Mesh->None
			,PlotLegends->{"mean angle = "<>ToString[DecimalForm[#,{3,1}]]<>"\[ThinSpace]\[Degree]"}
			,ScalingFunctions->{"Log10","Log10","Linear"}
			,Ticks->{trTicks,ttTicks,{1,5,10}}
		]&@meanAngleInterpolation;
		
	,
		(** fix for width outside of bounds issue (width capping) **)
		Print["Warning (i="<>ToString[i]<>"): Experimental width ("<>ToString[widthEx]<>") is outside the interval of numerically found widths "<>ToString[MinMax[cubeData[[All,4]]]]<>"! Setting parameter point to interval bounds now."];
		parameterPoint=If[widthEx>Max[cubeData[[All,4]]],parameterPointsStartEnd[[1]],parameterPointsStartEnd[[-1]]];
		width=widthInterpolationLogLogLin@@{Log10[#[[1]]],Log10[#[[2]]],#[[3]]}&@parameterPoint;
		meanAngle=angleInterpolationLogLogLin@@{Log10[#[[1]]],Log10[#[[2]]],#[[3]]}&@parameterPoint;
	];

	(** output **)
	{width,meanAngle,parameterPoint,plotWidth,plotAngle}

,{i,Length[widthsEx]}]\[Transpose]

]


(* ::Subsection::Closed:: *)
(*lineLinear*)


(* ::Text:: *)
(*linear interpolation between points start and end as a function of parameter t*)


lineLinear[t_,start_,end_]:=(1-t)*start + t*end;


(* ::Subsection::Closed:: *)
(*lineLinearLinearLinear*)


(* ::Text:: *)
(*linear change in linear space*)
(**)
(*inputs in : linear/linear/linear space*)
(*output : log/log/linear space*)


lineLinearLinearLinear[t_,parameterPointsStartEnd_]:={Log10[#[[1]]],Log10[#[[2]]],#[[3]]} &@ ( (1-t) parameterPointsStartEnd[[1]]+t parameterPointsStartEnd[[-1]] );


(* ::Subsection::Closed:: *)
(*lineLogLogLinear*)


(* ::Text:: *)
(*linear change in log-log-linear space*)
(**)
(*inputs in: linear/linear/linear space*)
(*output: log/log/linear space*)


lineLogLogLinear[t_,parameterPointsStartEnd_]:=((1-t){Log10[#[[1]]],Log10[#[[2]]],#[[3]]}&@ parameterPointsStartEnd[[1]]+t {Log10[#[[1]]],Log10[#[[2]]],#[[3]]}&@parameterPointsStartEnd[[-1]]);


(* ::Subsection:: *)
(*loadSimulationScanData: load parameter cube data for run time, tumble time, speed, intensity amplitude, width and mean angle*)


(* ::Input:: *)
(*dataScan=loadSimulationScanData[p];*)


(* ::Subsubsection:: *)
(*code*)


loadSimulationScanData[p_]:=Module[{simFileNames,simFilesSorted,datasets,extract,dataScan},
	simFileNames=FileNames[FileNameJoin[{p["path","parameterSweepNumerical"],"simAnalysisData_rDroplet_*.h5"}]];
	simFilesSorted=simFileNames[[Ordering[ToExpression[StringSplit[FileBaseName[#],"_"][[-1]]&/@simFileNames]]]];
	datasets=Import[#,"/Dataset1"]&/@simFilesSorted;
	
	(** output **)
	Flatten[Table[
		extract=GeneralUtilities`AssociationTranspose[dataset][#]&/@{"tRun","tTumble","u0Dim","iWidth","meanTheta","iminmax"};
		Transpose[({#[[1]],#[[2]],10^6#[[3]],#[[4]],180/\[Pi]#[[5]],#[[6]]})&@extract]
	,{dataset,datasets}],1]
]


(* ::Section:: *)
(*Plotting*)


(* ::Input:: *)
(*fig4[p];*)


(* ::Chapter::Closed:: *)
(*Figure S1: Liquid crystal optics*)


(* ::Section:: *)
(*Functions*)


(* ::Subsection:: *)
(*figS1*)


(* ::Input:: *)
(*figS1[p];*)


(* ::Subsubsection::Closed:: *)
(*code*)


figS1[p_]:=Module[{pthOutputFigS1,solution,pthIntensityImages,pthIntensityLookupMap},
	{pthOutputFigS1,solution,pthIntensityImages,pthIntensityLookupMap}=figS1PrepareData[p];
	figS1a[solution,pthOutputFigS1,p];
	figS1b[solution,pthOutputFigS1,p];
	figS1c[solution,pthOutputFigS1,p];
	(*figS1d[pthOutputFigS1,pthIntensityImages,3,3];*)
	figS1d[pthOutputFigS1,pthIntensityImages];
	figS1e[pthOutputFigS1,pthIntensityLookupMap];
	figS1f[pthOutputFigS1,pthIntensityLookupMap];
]


(* ::Subsection::Closed:: *)
(*figS1a: 3d droplet*)


(* ::Input:: *)
(*figS1a[solution,pthOutputFigS1,p];*)


(* ::Subsubsection::Closed:: *)
(*test: liquid crystal alignment*)


(* ::Input:: *)
(*Show[*)
(*	Graphics3D[{Opacity[0.2],Sphere[]}],*)
(*	renderLCs[solution,{0.3,0.5}]*)
(*]*)


(* ::Subsubsection::Closed:: *)
(*test: liquid crystal visualization*)


(* ::Input:: *)
(*liquidCrystals=Graphics3D[{Darker[Red],Evaluate[Cylinder[#,0.02]&/@vectors]}];*)
(*Show[hemisphereTop,liquidCrystals,hemisphereBottom,PlotRange->All,Lighting->"Neutral"]*)


(* ::Subsubsection::Closed:: *)
(*test: biphase droplet visualization *)


(* ::Input:: *)
(*hemisphereTop=ParametricPlot3D[{Cos[u] Sin[v],Sin[u] Sin[v],Cos[v]},{u,0,2\[Pi]},{v,0,\[Pi]/2},Mesh->None,Boxed->False,Axes->None,PlotPoints->50,PlotStyle->Directive[Darker[Red],Opacity[0.3]],Lighting->(*"Accent"*){White}];*)
(*hemisphereBottom=ParametricPlot3D[{Cos[u] Sin[v],Sin[u] Sin[v],Cos[v]},{u,0,2\[Pi]},{v,\[Pi]/2,\[Pi]},Mesh->None,Boxed->False,Axes->None,PlotPoints->50,PlotStyle->Directive[White,Opacity[0.2]],Lighting->"Accent"(*{White}*)];*)
(*Show[hemisphereTop,hemisphereBottom,PlotRange->All]*)


(* ::Subsubsection:: *)
(*code*)


figS1a[solution_,pthOutputFigS1_,p_]:=Module[{rs,n\[Theta]s,n\[Phi]s,r,n\[Theta],n\[Phi],r\[Theta]\[Phi]points,vectors,x,y,z,orientation,plotOps,plot,lcColor,hemisphereTop,liquidCrystals,hemisphereBottom},
	
	(** define r\[Theta]\[Phi]-sampling points **)
	rs={0.96,0.7,0.4,0.1};
	n\[Theta]s={7,5,3,2};
	n\[Phi]s={15,15,15,3};
	r\[Theta]\[Phi]points=Flatten[MapThread[
		({r,n\[Theta],n\[Phi]}={#1,#2,#3};
		r Flatten[
		Table[
			N[{#[[1]],#[[2]],Cos[\[Theta]]}&/@(Sin[\[Theta]]CirclePoints[Round[n\[Phi] Sin[\[Theta]]+5]])]
		,{\[Theta],Subdivide[0.15,0.95 \[Pi]/2,n\[Theta]-1]}]
		,1])&
	,{rs,n\[Theta]s,n\[Phi]s}],1];
	
	(** get liquid crystal directors **)
	Clear[u,v,w];
	vectors=Table[
		{x,y,z}=xyz;
		orientation=Normalize[{u[x,y,z],v[x,y,z],w[x,y,z]}/.solution];
		(xyz+# 0.05orientation)&/@{-1,1}
	,{xyz,r\[Theta]\[Phi]points}];
	
	(** plotting **)
	plotOps={Mesh->None,Boxed->False,Axes->None,PlotPoints->50};
	hemisphereTop=ParametricPlot3D[{Cos[u] Sin[v],Sin[u] Sin[v],Cos[v]},{u,0,2\[Pi]},{v,0,\[Pi]/2},Evaluate[plotOps],PlotStyle->Directive[Darker[Red],Opacity[0.3]],Lighting->(*"Accent"*){White}];
	hemisphereBottom=ParametricPlot3D[{Cos[u] Sin[v],Sin[u] Sin[v],Cos[v]},{u,0,2\[Pi]},{v,\[Pi]/2,\[Pi]},Evaluate[plotOps],PlotStyle->Directive[White,Opacity[0.2]],Lighting->"Accent"(*{White}*)];
	lcColor=RGBColor[Rational[71,85],Rational[42,85],Rational[42,85]];
	liquidCrystals=Graphics3D[{Opacity[0.8],lcColor,EdgeForm[None],Evaluate[Cylinder[#,0.02]&/@vectors]}];
	plot=Show[hemisphereTop,liquidCrystals,hemisphereBottom,PlotRange->All,Lighting->"Accent"(*"Neutral"*)];
	Print[plot];
	Export[FileNameJoin[{pthOutputFigS1,"fig_s1a_lc_droplet.png"}],plot];
]


(* ::Subsection::Closed:: *)
(*figS1b: xz-cut through volume *)


(* ::Input:: *)
(*figS1b[solution,pthOutputFigS1,p];*)


(* ::Subsubsection:: *)
(*code*)


figS1b[solution_,pthOutputFigS1_,p_]:=Module[{effectiveRadius,rMax,nxSamples,nzSamplesMax,xzpoints,zMin,zMax,dz,lcColor,liquidCrystals,hemisphereTop,hemisphereBottom,plot,x,y,z,yCut,orientation,vectors,plotOps},
	
	effectiveRadius=0.99;
	rMax=0.95;
	
	(** define xz-sampling points **)
	{nxSamples,nzSamplesMax}={10,5};
	xzpoints=Flatten[Table[
		zMax=Sqrt[effectiveRadius^2-x^2];
		zMin=0.05+1-effectiveRadius;
		dz=(zMax-zMin)/nzSamplesMax;
		If[dz<0.1,dz*=2];
		
		(** output **)
		If[zMax>zMin,
			Table[{x,z},{z,N@Range[zMin,zMax,dz]}]
			,Sequence@@{}
		]
		,{x,N[Subdivide[-rMax,rMax,nxSamples-1]]}
	],1];
	(*{nrSamples,n\[Phi]Samples}={10,10};
	points=Flatten[Table[
	Table[ {r Cos[\[Phi]],r Sin[\[Phi]]+0.1},{r,N@Subdivide[0.1,effectiveRadius,nrSamples-1]}]
	,{\[Phi],N[Subdivide[0,\[Pi],n\[Phi]Samples-1]]}],1];*)
	
	(** get liquid crystal directors **)
	yCut=0.0;
	Clear[u,v,w];
	vectors=Table[
		{x,y,z}={xz[[1]],yCut,xz[[2]]};
		orientation=Normalize[{u[x,y,z],v[x,y,z],w[x,y,z]}/.solution];
		(*orientation=directorAnalytical[{x,y,z}];*)
		(** output **)
		({x,y,z}+# 0.05orientation)&/@{-1,1}
		(*(p+# 0.05director[p])&/@{-1,1}*)
	,{xz,xzpoints}];
	
	(** plotting **)
	plotOps={Mesh->None,Boxed->False,Axes->None,PlotPoints->50};
	hemisphereTop=ParametricPlot3D[r{Cos[u],0,Sin[u]},{u,0,\[Pi]},{r,0,1},Evaluate[plotOps],PlotStyle->Directive[Darker[Red],Opacity[0.3]],Lighting->{White}];
	hemisphereBottom=ParametricPlot3D[r{Cos[u],0,Sin[u]},{u,\[Pi],2\[Pi]},{r,0,1},Evaluate[plotOps],PlotStyle->Directive[Lighter@Gray,Opacity[0.2]],Lighting->"Accent"];
	lcColor=RGBColor[Rational[71,85],Rational[42,85],Rational[42,85]];
	liquidCrystals=Graphics3D[{Opacity[0.9],lcColor,EdgeForm[None],Evaluate[Cylinder[#,0.02]&/@vectors]},PlotRange->{All,All,1.1{-1,1}},Boxed->False,ViewVertical->{0,0,1},ViewPoint->{0,\[Pi],0}];
	plot=Show[liquidCrystals,hemisphereTop,hemisphereBottom,PlotRange->1.05{{-1,1},{-1,1}}];
	Print[plot];
	Export[FileNameJoin[{pthOutputFigS1,"fig_s1b_lc_droplet_xz_cross_section.png"}],plot];
]


(* ::Subsection::Closed:: *)
(*figS1c: top view of curved interface *)


(* ::Input:: *)
(*figS1c[solution,pthOutputFigS1,p];*)


(* ::Subsubsection:: *)
(*code*)


figS1c[solution_,pthOutputFigS1_,p_]:=Module[{effectiveRadius,rMax,plotOps,nrSamples,n\[Phi]SamplesMin,n\[Phi]SamplesMax,r\[Phi]points,n\[Phi]Samples,vectors,x,y,z,orientation,liquidCrystals,hemisphereTop,plot,lcColor},
	
	effectiveRadius=0.99;
	rMax=0.95;
	
	(** define r\[Phi]-sampling points **)
	{nrSamples,n\[Phi]SamplesMin,n\[Phi]SamplesMax}={5,6,20};
	r\[Phi]points=Flatten[Table[
		n\[Phi]Samples=Round[n\[Phi]SamplesMin+(n\[Phi]SamplesMax-n\[Phi]SamplesMin)/(rMax-0.1)r];
		Table[r{Cos[\[Phi]],Sin[\[Phi]]},{\[Phi],N[Subdivide[0,2\[Pi],n\[Phi]Samples-1][[;;-2]]]}]
		,{r,N@Subdivide[0.1,rMax,nrSamples-1]}
	],1];
	
	(** get liquid crystal directors **)
	Clear[u,v,w];
	vectors=Table[
		{x,y,z}={r\[Phi][[1]],r\[Phi][[2]],Sqrt[effectiveRadius^2-r\[Phi][[1]]^2-r\[Phi][[2]]^2]};
		orientation=Normalize[{u[x,y,z],v[x,y,z],w[x,y,z]}/.solution];
		(*orientation=directorAnalytical[{x,y,z}];*)
		({x,y,0}+# 0.05orientation)&/@{-1,1}
	,{r\[Phi],r\[Phi]points}];

	
	(** plotting **)
	plotOps={Mesh->None,Boxed->False,Axes->None,PlotPoints->50};
	hemisphereTop=ParametricPlot3D[r{Cos[u],Sin[u],0}+{0,0,0.1},{u,0,2\[Pi]},{r,0,1},Evaluate[plotOps],PlotStyle->Directive[Darker[Red],Opacity[0.3]],Lighting->{White}];
	lcColor=RGBColor[Rational[71,85],Rational[42,85],Rational[42,85]];
	liquidCrystals=Graphics3D[{Opacity[0.9],lcColor,EdgeForm[None],Evaluate[Cylinder[#,0.02]&/@vectors]},PlotRange->{All,All,1.1{-1,1}}];
	plot=Show[liquidCrystals,hemisphereTop,PlotRange->1.05{{-1,1},{-1,1}},Boxed->False,ViewPoint->{0,0,\[Pi]},ViewVertical->{0,-1,0}(*,ViewVertical\[Rule]{0,0,1},ViewPoint\[Rule]{0,\[Pi],0}*)];
	Print[plot];
	Export[FileNameJoin[{pthOutputFigS1,"fig_s1c_lc_droplet_top_view.png"}],plot];
]


(* ::Subsection::Closed:: *)
(*figS1d: droplet intensity overview image*)


(* ::Input:: *)
(*figS1d[pthOutputFigS1,pthIntensityImages,3,3];*)


(* ::Text:: *)
(*Typically there are 40x20 images in theta and phi directions for Jones calculations and 16x8 for meep calculations*)


(* ::Subsubsection::Closed:: *)
(*history of paths to simulated intensity images*)


(* ::Input:: *)
(*pthIntensityImages=FileNameJoin[{NotebookDirectory[],"old output","matlab_output","sphericalwave_Hananah_7nm_test"}];*)
(*pthIntensityImages="C:\\Users\\Jan\\Dropbox (MIT)\\Bacteria Droplet Manuscript\\Code\\Droplet optics\\rgb avg for fig3\\RGB_averaged_images";*)
(*pthIntensityImages="C:\\Users\\Jan\\Dropbox (MIT)\\Bacteria Droplet Manuscript\\Code\\Droplet optics\\propagated_intensity_images";*)
(*pthIntensityImages=FileNameJoin[{p["path","dropletOpticsNumerical"],"intensity_images"}];*)


(* ::Subsubsection:: *)
(*code*)


Clear[figS1d];
figS1d[pthOutputFigS1_,pthImages_,dTheta_:1,dPhi_:1]:=Module[{pthNumericalIntensityImages,fileNames,filePrefix,intensityImageFiles,intensityImagesFlat,intensityImages,intensityOverviewImage,nTheta,nPhi},
	pthNumericalIntensityImages=FileNameJoin[{p["path","dropletOpticsNumerical"],"intensity_images"}];
	Print["Source for intensity images: ",pthNumericalIntensityImages];
	fileNames=FileNames[All,pthNumericalIntensityImages];
	filePrefix=StringSplit[FileBaseName[fileNames[[1]]],"_theta_"][[1]];
	{nTheta,nPhi}=Length/@DeleteDuplicates/@Transpose[StringSplit[FileBaseName[#],"_"][[{-3,-1}]]&/@fileNames];
	intensityImageFiles=FileNames[filePrefix<>"_theta_*_phi_*.png",pthImages];
	intensityImagesFlat=ImageReflect[ImageRotate[Import[#],-\[Pi]/2]]&/@intensityImageFiles;
	intensityImages=Partition[intensityImagesFlat,nPhi];
	intensityOverviewImage=ImageAssemble[Reverse[intensityImages\[Transpose][[;;;;dTheta,;;;;dPhi]]],Spacings->20,Background->White];
	Print[intensityOverviewImage];
	Export[FileNameJoin[{pthOutputFigS1,"fig_s1d_intensityOverviewImage_RGB_"<>ToString[dTheta]<>".png"}],intensityOverviewImage];
]


(* ::Subsection::Closed:: *)
(*figS1e: continuous mean intensity map*)


(* ::Input:: *)
(*figS1e[pthOutputFigS1,pthIntensityLookupMap];*)


(* ::Subsubsection:: *)
(*code*)


figS1e[pthOutputFigS1_,pthIntensityLookupMap_]:=Module[{spatiallyAveragedIntensities,plot},
	spatiallyAveragedIntensities=Import[FileNameJoin[{pthIntensityLookupMap,"normIntegratedIntensity.mat"}]][[1]];
	plot=ListDensityPlot[spatiallyAveragedIntensities\[Transpose]
		(*DataRange->{{0,90},{0,45}},*)
		,AspectRatio->1/2
		,ColorFunction->GrayLevel
		,PlotRangePadding->None
		,FrameStyle->Directive[Black,20,AbsoluteThickness[2]]
		,FrameLabel->{
			"Polar angle \!\(\*StyleBox[\"\[Theta]\",\nFontSlant->\"Italic\"]\) (\[Degree])",
			"Azimuthal angle \!\(\*StyleBox[\"\[Phi]\",\nFontSlant->\"Italic\"]\) (\[Degree])"
		}
		,FrameTicks->None
		,ImageSize->{Automatic,400}
	];
	Print[plot];
	Export[FileNameJoin[{pthOutputFigS1,"fig_s1e.png"}],plot];
]


(* ::Subsection::Closed:: *)
(*figS1f: continuous mean intensity map*)


(* ::Input:: *)
(*figS1f[pthOutputFigS1,pthIntensityLookupMap];*)


(* ::Subsubsection:: *)
(*code*)


figS1f[pthOutputFigS1_,pthIntensityLookupMap_]:=Module[{thetaPhiIntensityInterpolation,data,plot,plotArrowsQ},
	thetaPhiIntensityInterpolation=updateIntensityLookupMap[pthIntensityLookupMap];
	
	plotArrowsQ=False;
	
	(** \[Theta]=radial coordinate, \[Phi]=angle coordinate **)
	data=Flatten[Table[{\[Theta] Cos[\[Phi]],\[Theta] Sin[\[Phi]],thetaPhiIntensityInterpolation[-Abs[\[Theta]-\[Pi]/2.0]+\[Pi]/2.0,-Abs[Mod[\[Phi],\[Pi]/2.0]-\[Pi]/4.0]+\[Pi]/4.0]},{\[Theta],Subdivide[0,\[Pi]/2,100-1]},{\[Phi],Most[Subdivide[0,2\[Pi],100-1]]}],1];
	plot=ListDensityPlot[data
		,PlotRange->All
		,ColorFunction->GrayLevel
		,If[plotArrowsQ,
			Epilog->{Darker@Red,
				Arrow[{{0,0},{\[Pi]/2,0}}],Text[Style["\[Theta]",18],Scaled[{0.92,0.5}]],
				Arrow@ResourceFunction["SplineCircle"][{0,0},1.1\[Pi]/2,{1,0},{5Degree,90 Degree}],Text[Style["\[Phi]",18],Scaled[{0.83,0.83}]]
			}
			,Sequence@@{}
		]
		,Frame->False
		,PlotRangePadding->Scaled[0.1]
	];
	Print[plot];
	Export[FileNameJoin[{pthOutputFigS1,"fig_s1f.png"}],plot];
]


(* ::Subsection::Closed:: *)
(*figS1PrepareData*)


(* ::Input:: *)
(*{pthOutputFigS1,solution,pthIntensityImages,pthIntensityLookupMap}=figS1PrepareData[p];*)


(* ::Subsubsection:: *)
(*code*)


Clear[figS1PrepareData];
figS1PrepareData[p_]:=Module[{pthOutputFigS1,liquidCrystalSolution,pthIntensityImages,pthIntensityLookupMap},
	
	(** output directory **)
	pthOutputFigS1=FileNameJoin[{p["path","figs"],"figure S1"}];
	If[!DirectoryQ@#,CreateDirectory@#]&@pthOutputFigS1;
	
	(** solution for liquid crystal alignment **)
	liquidCrystalSolution=solveHemisphereLC[];
	
	(** path to simulated intensity images **)
	pthIntensityImages=FileNameJoin[{p["path","dropletOpticsNumerical"],"intensity_images"}];
	pthIntensityLookupMap=FileNameJoin[{p["path","dropletOpticsNumerical"],"lookup_map"}];
	
	(** output **)
	{pthOutputFigS1,liquidCrystalSolution,pthIntensityImages,pthIntensityLookupMap}
]


(* ::Section:: *)
(*Plotting*)


(* ::Input:: *)
(*figS1[p];*)


(* ::Chapter::Closed:: *)
(*Figure S2: Analysis pipeline*)


(* ::Subsection:: *)
(*figS2Prepare*)


(* ::Subsection::Closed:: *)
(*Figure S2c*)


(* ::Input:: *)
(*filterLists[[1]][[5]]*)


(* ::Input:: *)
(*pthExData=FileNameJoin[{ParentDirectory[NotebookDirectory[]],"tracker","output_2maxLength_frames"}];*)
(*filterQ=True;*)
(*badLists={(*{1,2,11,26,28,31,44,51,70}~Join~*){1}~Join~{3,7,9,10,16,20,23,27,29,41,43,44,46,49,50,54,55,59,64},{4,5,17,18,21,22,23,33,36,38,47,49},{1,2,3}~Join~{4,6,7,9,15,16,21,24,36,39,40,42,43,45,59,62,66,71,73},{3,6,7,11,13,14,15,24,30,31,32,33,40,43,44,51,57,65,77,79,80,82,85,92,96,98,103,105,108,114,117,118,124},*)
(*{5,13,21,32,33,34,52,55},{}};*)
(*numberDroplets={71,50,83,124,62,14};  (** number of droplets **)*)
(*filterLists=MapThread[Complement[Range[#2],#1]&,{badLists,numberDroplets}];  (** good droplets **)*)
(**)
(*allIntensities=Table[*)
(*intensityParts=Import[FileNameJoin[{pthExData,"all_intensityParts_"<>ToString[dayIndex]<>".dat"}]];*)
(*If[filterQ,intensityParts=intensityParts[[filterLists[[dayIndex]]]];];*)
(**)
(*(** output **)*)
(*intensityParts*)
(*,{dayIndex,Range[nDays]}];*)


(* ::Input:: *)
(*dt=1.0/30.0;*)
(*Dimensions/@({Range[Length[#]]dt,#}\[Transpose]&/@allIntensities[[1]])*)


(* ::Input:: *)
(*colorHues[col_]:=Flatten[{Lighter[col,#]&/@Reverse@Subdivide[0.1,0.7,10],col,Darker[col,#]&/@Subdivide[0.1,0.7,10]}]*)
(*colorHues[colorYoung]*)


(* ::Input:: *)
(*Length[allIntensities[[1]]]*)
(*iTicks=Table[{i,DecimalForm[i,{2,1}]},{i,0,0.5,0.1}];*)
(*iTicks2={#[[1]],""}&/@iTicks;*)
(*tTicks=Table[t,{t,0,300,50}];*)
(*tTicks2={#,""}&/@tTicks;*)
(*plot=ListLinePlot[{Range[Length[#]]dt,#}\[Transpose]&/@allIntensities[[1]],PlotRange->{{0,All},{0,0.3}},Frame->True,FrameStyle->Directive[Black,13,AbsoluteThickness[1.5]],FrameLabel->{"Time \!\(\**)
(*StyleBox[\"t\",\nFontSlant->\"Italic\"]\) (s)","Integrated intensity \!\(\**)
(*StyleBox[\"I\",\nFontSlant->\"Italic\"]\)"},PlotRangePadding->None,FrameTicks->{{iTicks,iTicks2},{tTicks,tTicks2}},AspectRatio->0.56,(*PlotStyle->Directive[Opacity[0.2],colorYoung]*)PlotStyle->Reverse[colorHues[colorYoung]]]*)
(*Export[FileNameJoin[{"C:\\Users\\Jan\\Dropbox (MIT)\\Bacteria Droplet Manuscript\\Figures\\Figures current\\figure_S2_pipeline","fig_s2c.png"}],plot]*)


(* ::Input:: *)
(*Manipulate[*)
(*ListLinePlot[{Range[Length[#]]dt,#}\[Transpose]&@allIntensities[[1,i]],PlotRange->{{0,300},{0,0.3}},Frame->True,FrameStyle->Directive[Black,13,AbsoluteThickness[1.5]],FrameLabel->{"Time \!\(\**)
(*StyleBox[\"t\",\nFontSlant->\"Italic\"]\) (s)","Integrated intensity \!\(\**)
(*StyleBox[\"I\",\nFontSlant->\"Italic\"]\)"},PlotRangePadding->None,FrameTicks->{{iTicks,iTicks2},{tTicks,tTicks2}},AspectRatio->0.56,PlotStyle->Directive[Opacity[0.5],colorYoung](*PlotStyle->Reverse[colorHues[colorYoung]]*)]*)
(*,{i,1,Length[allIntensities[[1]]],1}]*)


(* ::Subsection::Closed:: *)
(*Figure S2d*)


(* ::Input:: *)
(*windowingQ=True;*)
(*fMax=1.0;*)
(**)
(*movieLengths={8344, 7252,8046,9225,5244,2251};*)
(*padFactor=2;*)
(*nPaddedSamples=padFactor Max[movieLengths]; *)
(*filterQ=True;*)
(*badLists={(*{1,2,11,26,28,31,44,51,70}~Join~*){1}~Join~{3,7,9,10,16,20,23,27,29,41,43,44,46,49,50,54,55,59,64},*)
(*{4,5,17,18,21,22,23,33,36,38,47,49},*)
(*{4,6,7,9,15,16,21,24,36,39,40,42,43,45,59,62,66,71,73},{3,6,7,11,13,14,15,24,30,31,32,33,40,43,44,51,57,65,77,79,80,82,85,92,96,98,103,105,108,114,117,118,124},*)
(*{5,13,21,32,33,34,52,55},*)
(*{}*)
(*};*)
(*numberDroplets={71,50,83,124,62,14};  (** number of droplets **)*)
(*filterLists=MapThread[Complement[Range[#2],#1]&,{badLists,numberDroplets}];  (** good droplets **)*)
(**)
(*(** calculate mean spectrum for each day **)*)
(*Print[AbsoluteTiming[*)
(*days={1};*)
(*{allIntensities,meanSpectra,stdSpectra,allSpectra}=Table[*)
(*intensityParts=Import[FileNameJoin[{pthData,"all_intensityParts_"<>ToString[dayIndex]<>".dat"}]];*)
(*If[filterQ,intensityParts=intensityParts[[filterLists[[dayIndex]]]];];*)
(**)
(*(** get mean/std spectrum **)*)
(*times=Range[0,Length[intensityParts[[1]]]-1]dt;*)
(*paddedMemory=ConstantArray[0.0,nPaddedSamples];*)
(*powerSpectralDensities=calculateTemporalFreqSpectrum[dt,#,windowingQ,nPaddedSamples,paddedMemory][[{1,2}]]\[Transpose]&/@intensityParts;*)
(*freqs=powerSpectralDensities[[1,All,1]];*)
(*spectra=powerSpectralDensities[[All,All,2]];*)
(*meanSpectrum=Mean[spectra];*)
(*stdSpectrum=StandardDeviation[spectra];*)
(**)
(*(** output **)*)
(*{intensityParts,meanSpectrum,stdSpectrum,spectra}*)
(*,{dayIndex,days}]\[Transpose];*)
(*][[1]]];*)


(* ::Input:: *)
(*fTicks=Table[If[Mod[t,0.5]==0,{t,DecimalForm[t,{2,1}],{0.03,0}},{t,"",{0.015,0}}],{t,0,1,0.1}];*)
(*fTicks2={#[[1]],"",#[[3]]}&/@fTicks;*)
(*{plotStdQ,plotIndiQ}={True,False};*)
(*pTicks=Table[{t,DecimalForm[t,{2,1}]},{t,0,1,0.5}];*)
(*pTicks2={#[[1]],""}&/@pTicks;*)
(*yMin=-0.02;*)
(*yMax=4.5;*)
(*fwhmPoints={{0.15,0.5},{0.07,0.5}};*)
(*space=-8;*)
(*xLabel=Framed["Frequency \!\(\**)
(*StyleBox[\"f\",\nFontSlant->\"Italic\"]\) (Hz)",FrameStyle->None,ContentPadding->False,FrameMargins->{{0,0},{0,space}}];*)
(*ops={PlotRange->{{-0.01,1.01}fMax,{yMin,yMax}},Frame->True,AspectRatio->1,FrameLabel->{xLabel,"Intensity spectrum"},FrameStyle->Directive[Black,fontSize,AbsoluteThickness[lineThickness],"Arial"],FrameTicks->{{pTicks,pTicks2},{fTicks,fTicks2}}};*)
(*dayIndices={1,6};*)
(**)
(*plotIndex=1;*)
(*color=If[plotIndex==1,colorYoung,colorOld];*)
(*meanSpectrum=meanSpectra[[dayIndices[[plotIndex]]]];*)
(*(** fit spectrum **)*)
(*fMax=1.0;*)
(*fMaxIndex=FirstPosition[freqs,_?(#>=fMax&)][[1]];*)
(*{indexStart,indexEnd}={1,fMaxIndex};*)
(*iFit=NonlinearModelFit[{freqs,#/Max[#]}\[Transpose][[indexStart;;indexEnd]],{aa/(1+(\[Omega]/bb)^2),aa>0,bb>0},{{aa,1},{bb,0.1}},\[Omega],MaxIterations->1000]&@meanSpectrum;   (** Lorentz with peak at 0 Hz **)*)
(*{iAmplitude,iWidth}=iFit["BestFitParameters"][[All,2]];*)
(*fwhmPoint={iWidth,0.5};*)
(**)
(*psdMax=iAmplitude Max[meanSpectrum];*)
(*meanPSD={freqs,#/psdMax}\[Transpose]&@meanSpectrum;*)


(* ::Input:: *)
(*pMeanAndIndividuals=Show[*)
(*ListLinePlot[{freqs,#/psdMax}\[Transpose],PlotRange->{{-0.01,1.01}fMax,{yMin,yMax}},FrameStyle->Directive[(*FontFamily->"Arial",*)Black,10,AbsoluteThickness[1.2]],AspectRatio->0.75,Evaluate[ops],PlotStyle->Directive[Opacity[0.1],colorYoung]]&/@allSpectra[[dayIndices[[plotIndex]]]],*)
(*ListLinePlot[meanPSD,PlotStyle->Directive[AbsoluteThickness[3],colorYoung],Evaluate[ops]](*Epilog->{Red,AbsolutePointSize[7.5],Point[#],AbsoluteThickness[lineThickness],Dashed,HalfLine[{#,{0,0.5}}],HalfLine[{#,{#\[LeftDoubleBracket]1\[RightDoubleBracket],0}}]}&@fwhmPoints\[LeftDoubleBracket]plotIndex\[RightDoubleBracket],*)*)
(*,FrameTicks->{{None,None},{fTicks,fTicks2}}*)
(*,ImageSize->{Automatic,45 mm}*)
(*(*,FrameStyle->Directive[(*FontFamily->"Arial",*)Black,10,AbsoluteThickness[1.2]]*)*)
(*]*)
(*Export[FileNameJoin[{"C:\\Users\\Jan\\Dropbox (MIT)\\Bacteria Droplet Manuscript\\Figures\\Figures current\\SI figures\\SI_figure_2_pipeline","fig_s2d.png"}],pMeanAndIndividuals]*)


(* ::Chapter:: *)
(*Figure S4: Mapping between widths and angles*)


(* ::Section:: *)
(*Functions*)


(* ::Subsection:: *)
(*figS4*)


(* ::Input:: *)
(*figS4[p];*)


(* ::Subsubsection:: *)
(*debug*)


	uMax=10;
	{widthsEx,dataScan,plotOpsFigS4,pthOutFigS4,meshPointsMinWidth,meshPointsMaxWidth,{wMin,wMax},filteredMappingsDec,meanMap,stdMap}=figS4Prepare[uMax,p];


(* ::Subsubsection:: *)
(*code*)


figS4[p_]:=Module[{(*uMax,widthsEx,dataScan,plotOpsFigS4,pthOutFigS4,meshPointsMinWidth,meshPointsMaxWidth,wMin,wMax,filteredMappingsDec,meanMap,stdMap*)},
	uMax=10;
	{widthsEx,dataScan,plotOpsFigS4,pthOutFigS4,meshPointsMinWidth,meshPointsMaxWidth,{wMin,wMax},filteredMappingsDec,meanMap,stdMap}=figS4Prepare[uMax,p];
	
	figS4a[dataScan,widthsEx,uMax,plotOpsFigS4,pthOutFigS4,p];
	figS4b[dataScan,widthsEx,uMax,plotOpsFigS4,pthOutFigS4,p];
	figS4c[dataScan,filteredMappingsDec,wMin,wMax,meshPointsMaxWidth,meshPointsMinWidth,uMax,plotOpsFigS4,pthOutFigS4,p];
	figS4d[filteredMappingsDec,meanMap,stdMap,pthOutFigS4,p];
]


(* ::Subsection:: *)
(*figS4a: inspect relevant width iso-surfaces*)


(* ::Text:: *)
(*plot width iso-surfaces in parameter space spanned by run-time, tumble-time and velocity*)


(* ::Input:: *)
(*figS4a[dataScan,widthsEx,uMax,plotOpsFigS4,pthOutFigS4,p];*)


(* ::Subsubsection::Closed:: *)
(*test: visually inspect iso-surfaces of width, angle and intensity fluctuations amplitude*)


(* ::Input:: *)
(*ops={ScalingFunctions->{"Log10","Log10"},ContourStyle->Opacity[0.6],Mesh->None,PlotLegends->Automatic,PlotRange->{Log10@{0.009,9.`},Log10@{0.01,1},{1,12}},AxesLabel->{"\!\(\*SubscriptBox[\(t\), \(Run\)]\) (s)","\!\(\*SubscriptBox[\(t\), \(Tumble\)]\) (s)",Rotate["u (\[Mu]m/s)",\[Pi]/2]},AxesStyle->Black};*)
(*ListContourPlot3D[dataScan[[All,{1,2,3,4}]],Evaluate[ops],PlotLabel->Style["Width (Hz)",Black],Contours->{0.08,0.161,0.24,0.25(*0.28*)}]*)
(*ListContourPlot3D[dataScan[[All,{1,2,3,5}]],Evaluate[ops],PlotLabel->Style["Mean angle (\[Degree])",Black],Contours->Range[5,25,10.]]*)
(*ListContourPlot3D[dataScan[[All,{1,2,3,6}]],Evaluate[ops],PlotLabel->"\!\(\*SubscriptBox[\(I\), \(max\)]\)-\!\(\*SubscriptBox[\(I\), \(min\)]\)"]*)


(* ::Input:: *)
(*ListContourPlot3D[dataScan[[All,{1,2,3,4}]],Evaluate[ops],PlotLabel->"Width (Hz)",Contours->{0.08,0.24}]*)


(* ::Subsubsection:: *)
(*code*)


figS4a[dataScan_,widthsEx_,uMax_,plotOpsFigS4_,pthOutFigS4_,p_]:=Module[{plot,meanAngleColors},
	
	plot=ListContourPlot3D[dataScan[[All,{1,2,3,4}]]
		,Evaluate[plotOpsFigS4]
		,Contours->(Round[#,0.001]&/@widthsEx)
		,ContourStyle->(Directive[Opacity[0.4],#]&/@p["color","ages"])
		,Lighting->"Accent"
		,PlotLabel->Style["Width (Hz)",Black]
		,PlotLegends->Automatic
	];
	Print[plot];
	Export[FileNameJoin[{pthOutFigS4,"fig_s4a.png"}],plot,ImageResolution->300];
]


(* ::Subsection::Closed:: *)
(*figS4b: inspect relevant mean angle iso-surfaces*)


(* ::Input:: *)
(*figS4b[dataScan,widthsEx,uMax,plotOpsFigS4,pthOutFigS4,p];*)


(* ::Subsubsection:: *)
(*code*)


figS4b[dataScan_,widthsEx_,uMax_,plotOpsFigS4_,pthOutFigS4_,p_]:=Module[{plot,meanAngleColors},
		
	meanAngleColors=ColorData["SunsetColors"][#]&/@Subdivide[0,0.7,4-1];
	plot=ListContourPlot3D[dataScan[[All,{1,2,3,5}]]
		,Evaluate[plotOpsFigS4]
		,Contours->Range[5,35,10]
		,ContourStyle->(Directive[Opacity[0.4],#]&/@meanAngleColors)
		,Lighting->{White}
		,PlotLabel->Style["Mean angle \[LeftAngleBracket]\[Theta]\[RightAngleBracket] (\[Degree])",Black]
		,PlotLegends->Automatic
	];
	Print[plot];
	Export[FileNameJoin[{pthOutFigS4,"fig_s4b.png"}],plot,ImageResolution->300]
]


(* ::Subsection:: *)
(*figS4c: evolution lines and surfaces*)


(* ::Input:: *)
(*figS4c[dataScan,filteredMappingsDec,wMin,wMax,meshPointsMaxWidth,meshPointsMinWidth,uMax,plotOpsFigS4,pthOutFigS4,p];*)


(* ::Subsubsection:: *)
(*code*)


figS4c[dataScan_,filteredMappingsDec_,wMin_,wMax_,meshPointsMaxWidth_,meshPointsMinWidth_,uMax_,plotOpsFigS4_,pthOutFigS4_,p_]:=Module[{widthInterpolationLinLinLin,mappings,grayColors,plot},
	
	widthInterpolationLinLinLin=getScanDataInterpolantsLinLinLin[dataScan][[1]];
	mappings=filteredMappingsDec[[;;;;50]];
	(*grayColors=GrayLevel/@Subdivide[0.2,0.7,Length[mappings]-1];*)
	grayColors=Hue[0,0.1,#]&/@Subdivide[0.3,0.9,Length[mappings]-1];
	
	plot=Show[
		(** surfaces **)
		ContourPlot3D[widthInterpolationLinLinLin[x,y,z],{x,0.1,10},{y,0.01,0.1},{z,1,10}
			,PlotLabel->Style["Trajectories",Black]
			,PlotRange->{{0.01,10},{0.01,0.1},{1,uMax}}
			,PlotRangePadding->Scaled[0.02]
			,Evaluate[plotOpsFigS4]
			,BoundaryStyle->None
			,ContourStyle->Evaluate[Directive[#,Opacity[0.6]]&/@{p["color","old"],p["color","young"]}]
			,Contours->{wMin,wMax}
			,PlotLegends->Table["width = "<>ToString[DecimalForm[w,{2,2}]]<>" Hz",{w,{wMin,wMax}}]
		]
		(** sample points **)
		,ListPointPlot3D[{meshPointsMaxWidth,meshPointsMinWidth}
			,BoxRatios->{1,1,1}
			,PlotStyle->{p["color","young"],p["color","old"]}
		]
		(** connection lines **)
		,Graphics3D[{Opacity[0.5],AbsoluteThickness[0.7],MapThread[{#1,Line[#2[[{1,2}]]]}&,{grayColors,mappings}]}]
	];
	Print[plot];
	Export[FileNameJoin[{pthOutFigS4,"fig_s4c.png"}],plot,ImageResolution->300];
]


(* ::Subsection::Closed:: *)
(*figS4d: individual mappings*)


(* ::Input:: *)
(*figS4d[filteredMappingsDec,meanMap,stdMap,pthOutFigS4,p];*)


(* ::Subsubsection:: *)
(*code*)


figS4d[filteredMappingsDec_,meanMap_,stdMap_,pthOutFigS4_,p_]:=Module[{wTicks,wTicks2,\[Theta]Ticks,\[Theta]Ticks2,mappings,grayColors,stdUpMap,stdDownMap,plot},
	wTicks=Range[0.05,0.25,0.05];
	wTicks2={#,""}&/@wTicks;
	\[Theta]Ticks=Range[0,40,10];
	\[Theta]Ticks2={#,""}&/@\[Theta]Ticks;
	
	mappings=filteredMappingsDec[[All,3]][[;;;;50]];
	(*grayColors=GrayLevel/@Subdivide[0.2,0.7,Length[mappings]-1];*)
	grayColors=Hue[0,0.1,#]&/@Subdivide[0.3,0.9,Length[mappings]-1];
	
	stdUpMap={meanMap[[All,1]],meanMap[[All,2]]+stdMap[[All,2]]}\[Transpose];
	stdDownMap={meanMap[[All,1]],meanMap[[All,2]]-stdMap[[All,2]]}\[Transpose];
	
	
	plot=Show[
		ListLinePlot[mappings
			,Evaluate[p["plot","ops"]]
			,ImagePadding->{{60,11},{45,10}}
			,FrameTicks->{{\[Theta]Ticks,\[Theta]Ticks2},{wTicks,wTicks2}}
			,PlotRangeClipping->False
			,PlotStyle->(Directive[#,Opacity[0.5]]&/@grayColors)
			,AspectRatio->1
			,Epilog->{
				(** frame labels **)
				Text[Style[p["plot","widthLabel"],p["plot","fontSize"]],Scaled[{0.5,-0.12}]]
				,Rotate[Text[Style["Angle \[LeftAngleBracket]\[Theta]\[RightAngleBracket] (\[Degree])",p["plot","fontSize"]],Scaled[{-0.13,0.5}]],90\[Degree]]
			}
		]
		(** mean mapping **)
		,ListLinePlot[Sort/@{meanMap,stdUpMap,stdDownMap}
			,PlotStyle->(Directive[Black,(*Dotted*)(*Dashing[{0.02,0.01}(*{0.1,0.05}*)],*)AbsoluteThickness[# p["plot","lineThickness"]]]&/@{3,1.5,1.5})
			,Filling->{2->{3}}
		]
		
	];
	Print[plot];
	Export[FileNameJoin[{pthOutFigS4,"fig_s4d.png"}],plot,ImageResolution->300];
]


(* ::Subsection:: *)
(*figS4Prepare*)


(* ::Subsubsection:: *)
(*debug*)


(* ::Input:: *)
(*widthsEx=Flatten[Import[FileNameJoin[{p["path","widthsExperiment"],"experimental_widths.dat"}]]];*)
(*{dataScan,meshPointsMinWidth,meshPointsMaxWidth,{wMin,wMax},widthInterpolationLogLogLin,angleInterpolationLogLogLin}=getSampledWidthAngleContours[widthsEx,p];*)


(* ::Subsubsection:: *)
(*code*)


(* ::Input:: *)
(*{widthsEx,dataScan,plotOpsFigS4,pthOutFigS4,meshPointsMinWidth,meshPointsMaxWidth,{wMin,wMax},filteredMappingsDec}=figS4Prepare[uMax,p];*)


figS4Prepare[uMax_,p_]:=Module[{widthsEx,dataScan,plotOpsFigS5,trTicks,ttTicks,plotOpsFigS4,pthOutFigS4,meshPointsMinWidth,meshPointsMaxWidth,wMin,wMax,widthInterpolationLogLogLin,angleInterpolationLogLogLin,filteredMappingsDec,meanMap,stdMap},

	(** get previously calculated widths from experiment **)
	widthsEx=Flatten[Import[FileNameJoin[{p["path","widthsExperiment"],"experimental_widths.dat"}]]];
	
	(** plot options **)
	{trTicks,ttTicks}=setRunTumbleTicks[];
	plotOpsFigS4={
		AxesLabel->{
			"\!\(\*SubscriptBox[\(t\), \(Run\)]\) (s)"
			,""
			,Rotate["u (\[Mu]m/s)",\[Pi]/2]
		}
		,AxesStyle->Directive[Black,11]
		,BoundaryStyle->None
		,BoxStyle->Black
		,Epilog->{
			Text["\!\(\*SubscriptBox[\(t\), \(Tumble\)]\) (s)", Scaled[{0.16, 0.88}]]
		}
		,ImagePadding->{{30,0},{15,5}}
		,ImageSize->220
		,Mesh->None
		,PlotRange->{Log10@{0.01,10},Log10@{0.01,(*1*)0.1},{1,uMax}}
		,ScalingFunctions->{"Log10","Log10","Linear"}
		,Ticks->{trTicks,ttTicks,{1,5,10}}
	};
	
	pthOutFigS4=FileNameJoin[{p["path","figs"],"figure S4"}];
	If[!DirectoryQ@#,CreateDirectory@#]&@pthOutFigS4;
	
	(** get simulation scan data and more **)
	{dataScan,meshPointsMinWidth,meshPointsMaxWidth,{wMin,wMax},widthInterpolationLogLogLin,angleInterpolationLogLogLin}=getSampledWidthAngleContours[widthsEx,p];
	{filteredMappingsDec,meanMap,stdMap}=getWidthAngleMappings[meshPointsMinWidth, meshPointsMaxWidth,wMax,widthInterpolationLogLogLin,angleInterpolationLogLogLin][[;;3]];
	
	{widthsEx,dataScan,plotOpsFigS4,pthOutFigS4,meshPointsMinWidth,meshPointsMaxWidth,{wMin,wMax},filteredMappingsDec,meanMap,stdMap}
]


(* ::Section:: *)
(*Subfigures*)


(* ::Input:: *)
(*figS4[p];*)


(* ::Subsection::Closed:: *)
(*(old) intersections between evolution line and surfaces*)


(* ::Input:: *)
(*iSave=2;*)
(*showIntersectionPointsQ=True;*)
(**)
(*plot=Show[*)
(*plotWidthSave*)
(*,plotAngleSave*)
(*,Graphics3D[{*)
(*(** intersection points **)*)
(*	If[showIntersectionPointsQ,MapThread[{#2,Sphere[{Log10[#1[[1]]],Log10[#1[[2]]],#1[[3]]},{0.05,0.05,1.^3}]}&,{parameterPoints,p["color","ages"]}],Nothing],*)
(*(** end points **)*)
(*{Black,Sphere[{Log10[#[[1]]],Log10[#[[2]]],#[[3]]},{0.05,0.05,1.^3}]&/@(*parameterPoints[[{1,-1}]]*)parameterPointsStartEnd},*)
(*(** curve between first and last point **)*)
(*{Black,Line[(*{Log10[#\[LeftDoubleBracket]1\[RightDoubleBracket]],Log10[#\[LeftDoubleBracket]2\[RightDoubleBracket]],#\[LeftDoubleBracket]3\[RightDoubleBracket]}&/@*)(hypotheticalEvolution[#]&/@Subdivide[0.0,1,100])]}*)
(*(** selected point **)*)
(*(*,{Red,Sphere[{Log10[#\[LeftDoubleBracket]1\[RightDoubleBracket]],Log10[#\[LeftDoubleBracket]2\[RightDoubleBracket]],#\[LeftDoubleBracket]3\[RightDoubleBracket]},{0.05,0.05,1.^3}]&@parameterPoints[[iSave]]}*)*)
(*}]*)
(*,AxesLabel->{"\!\(\*SubscriptBox[\(t\), \(Run\)]\) (s)","\!\(\*SubscriptBox[\(t\), \(Tumble\)]\) (s)",Rotate["u (\[Mu]m/s)",\[Pi]/2]}*)
(*,AxesStyle->Directive[Black,11]*)
(*,BoxStyle->Directive[Black]*)
(*,BoundaryStyle->None*)
(*,ImageSize->220*)
(*,Lighting->{White}*)
(*,PlotRange->{Log10@{0.09,9.`},Log10@{0.01,1},{1,uMax}}*)
(*,ScalingFunctions->{"Log10","Log10","Linear"}*)
(*]*)
(*Export[FileNameJoin[{pthOutFigS5,"fig_s5d.png"}],plot]*)
