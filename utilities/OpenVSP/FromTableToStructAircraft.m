%% OPENVsp Aero Analysis Output into AIRCRAFT DATA STRUCTURE
% In this script we use the Open VSP output to obtain the lift and drag
% coefficient distribution along the main wing span, in correspondence to a
% resultant lift coefficient equal to one.

%integral coefficients
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.OpenVSP.CL.value    = res.CL;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.OpenVSP.CDo.value   = res.CDo;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.OpenVSP.CDi.value   = res.CDi;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.OpenVSP.CDtot.value = res.CDtot;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.OpenVSP.CS.value    = res.CS;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.OpenVSP.E.value     = res.E;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.OpenVSP.e.value     = res.e;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.OpenVSP.CFx.value   = res.CFx;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.OpenVSP.CFy.value   = res.CFy;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.OpenVSP.CFz.value   = res.CFz;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.OpenVSP.CMx.value   = res.CMx;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.OpenVSP.CMy.value   = res.CMy;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.OpenVSP.CMz.value   = res.CMz;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.OpenVSP.CMl.value   = res.CMl;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.OpenVSP.CMm.value   = res.CMm;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.OpenVSP.CMn.value   = res.CMn;

%loads
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.OpenVSP.Yavg.value    = res.Yavg;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.OpenVSP.Cl.value      = res.Cl;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.OpenVSP.Cd.value      = res.Cd;
Aircraft.Certification.Regulation.SubpartC.Flightloads.Balancingloads.OpenVSP.Cmy.value     = res.Cmy;



