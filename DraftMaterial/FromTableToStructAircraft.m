%% LIFT AND DRAG ALONG THE SPAN - CL = 1 
% In this script we use the Open VSP output to obtain the lift and drag
% coefficient distribution along the main wing span, in correspondence to a
% resultant lift coefficient equal to one.

res_dir = [Aircraft.Certification.Aircraft_Name.value '_results'];
cd(res_dir)
polar   = importdata([Aircraft.Certification.Aircraft_Name.value '_DegenGeom.polar']);
Aircraft.Certification.Regulation.SubpartC.Balancingloads.OpenVSP.CL    = polar.data(:,5);
Aircraft.Certification.Regulation.SubpartC.Balancingloads.OpenVSP.CDo   = polar.data(:,6);
Aircraft.Certification.Regulation.SubpartC.Balancingloads.OpenVSP.CDi   = polar.data(:,7);
Aircraft.Certification.Regulation.SubpartC.Balancingloads.OpenVSP.CDtot = polar.data(:,8);
Aircraft.Certification.Regulation.SubpartC.Balancingloads.OpenVSP.CS    = polar.data(:,9);
Aircraft.Certification.Regulation.SubpartC.Balancingloads.OpenVSP.E     = polar.data(:,10);
Aircraft.Certification.Regulation.SubpartC.Balancingloads.OpenVSP.e     = polar.data(:,11);
Aircraft.Certification.Regulation.SubpartC.Balancingloads.OpenVSP.CFx   = polar.data(:,12);
Aircraft.Certification.Regulation.SubpartC.Balancingloads.OpenVSP.CFy   = polar.data(:,13);
Aircraft.Certification.Regulation.SubpartC.Balancingloads.OpenVSP.CFz   = polar.data(:,14);
Aircraft.Certification.Regulation.SubpartC.Balancingloads.OpenVSP.CMx   = polar.data(:,15);
Aircraft.Certification.Regulation.SubpartC.Balancingloads.OpenVSP.CMy   = polar.data(:,16);
Aircraft.Certification.Regulation.SubpartC.Balancingloads.OpenVSP.CMz   = polar.data(:,17);
Aircraft.Certification.Regulation.SubpartC.Balancingloads.OpenVSP.CMl   = polar.data(:,18);
Aircraft.Certification.Regulation.SubpartC.Balancingloads.OpenVSP.CMm   = polar.data(:,19);
Aircraft.Certification.Regulation.SubpartC.Balancingloads.OpenVSP.CMn   = polar.data(:,20);