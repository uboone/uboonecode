from ROOT import TFile, TH1F, TCanvas, gDirectory, gROOT, gPad, TLegend, TPad, TLine, TH1D
import math

gROOT.ForceStyle()

### User-defined variables ###

# ROOT files with histograms
snfile = TFile("/home/jcrespo/MicroBooNE/SNAnalysis/19021/spectra/SNRun19021_960files3planes_SNMichelAna_hist_v5samxrd.root")
extunbfile = TFile("/home/jcrespo/MicroBooNE/SNAnalysis/EXTUNB/spectra/Run3_25kfiles3planes_SNMichelAna_hist_v5.root")
#extunbfile_emu = TFile("/home/jcrespo/MicroBooNE/SNAnalysis/Emulation/my_HitMultiplicity_IntegralCollaborationMeetingAllEvents.root")
extunbfile_emu = None
#extunbfile_emumerged = TFile("/home/jcrespo/MicroBooNE/SNAnalysis/Emulation/my_CorrectMergedHitsAlgorithm.root")
extunbfile_emumerged = None

# Name of histograms to retrieve
hlistname = [ "heSpectrum",
          "hgSpectrum",
          "htotSpectrum",
          "heHitMult",
          "hgHitMult",
          "htotHitMult",
          "heHitSpectrum",
          "hgHitSpectrum",
          "htotHitSpectrum",
          "heAngle",
          "heLength",
          "heLengthW",
          "heLengthT",
          "hgClusMult",
          "hgClusSpectrum",
          "hgClusHitMult",
          "hEventHitMult",
          "hEventHitSpectrum" ]

# Fraction of canvas for the bottom pad (make it 0 if you do not want one)
frontierpad = 0.3
#frontierpad = 0
logx = False
#logx = True
logy = False
#logy = True

# Event sizes in number of samples
snevtsize = 3200. # Only unique samples
extunbevtsize = 6400.
extunbevtsize_emu = 9595.

### End of user-defined variables ###

# Lists to hold histograms retrieved from files
snhlist= []
extunbhlist = []
extunbhdict_emu = {}
extunbhdict_emumerged = {}

# Normalization factor for the 3 TPC planes
# To do: change to single number (all planes had the same number of events)
normfactor = [0., 0., 0.]
normfactor_emu = 0.
normfactor_emumerged = 0.

# Retrieve objects from SN stream file
for rootdir in snfile.GetListOfKeys():
    snfile.cd( rootdir.GetName() )
    for hname in hlistname:
        htemp = TH1F()
        gDirectory.GetObject( hname, htemp )
        htemp.SetNameTitle( "sn_" + htemp.GetName() + "_" + rootdir.GetName()[-1], htemp.GetTitle() + " plane " + rootdir.GetName()[-1] )
        htemp.Sumw2()
        snhlist.append( htemp )
        # Get normalization factor numerator
        if htemp.GetName().startswith("sn_hEventHitMult"):
            # print "%s %i " % (htemp.GetName(), htemp.GetEntries())
            normfactor[int( rootdir.GetName()[-1] )] = htemp.GetEntries()*snevtsize
            snexposure = htemp.GetEntries()*snevtsize*0.5e-6 # in seconds
            print "SN stream exposure: %i events = %.2f min " % (htemp.GetEntries(), snexposure/60.)
            if normfactor_emu == 0.: # Do it only once (all planes had the same number of events)
                normfactor_emu = htemp.GetEntries()*snevtsize
            if normfactor_emumerged == 0.: # Do it only once (all planes had the same number of events)
                normfactor_emumerged = htemp.GetEntries()*snevtsize

# Retrieve objects from EXT unbiased file
for rootdir in extunbfile.GetListOfKeys():
    extunbfile.cd( rootdir.GetName() )
    for hname in hlistname:
        htemp = TH1F()
        gDirectory.GetObject( hname, htemp )
        # Get exposure (number of events) from hEventHitMult
        htemp.SetNameTitle( "extunb_" + htemp.GetName() + "_" + rootdir.GetName()[-1], htemp.GetTitle() + " plane " + rootdir.GetName()[-1] )
        # Get normalization factor denominator
        if htemp.GetName().startswith("extunb_hEventHitMult"):
            # print "%s %i " % (htemp.GetName(), 2.*htemp.GetEntries())
            normfactor[int( rootdir.GetName()[-1] )] *= 1./(htemp.GetEntries()*extunbevtsize)
            extunbexposure = htemp.GetEntries()*extunbevtsize*0.5e-6 # in seconds
            print "Trigger stream exposure: %i events = %.2f min " % (htemp.GetEntries(), extunbexposure/60.)
            # For hit multiplicity/event we have to correct for the different event size
            hnew = htemp.Clone()
            hnew.Reset()
            for ibin in xrange( htemp.GetNbinsX() ):
                # Trigger stream events are twice as long as SN stream ones, so they can have twice as many hits
                ibinnew = int( math.floor(ibin*snevtsize/extunbevtsize) )
                hnew.SetBinContent( ibinnew, hnew.GetBinContent(ibinnew) + htemp.GetBinContent(ibin) )
            htemp = hnew    
        htemp.Sumw2()
        extunbhlist.append( htemp )

# Retrieve objects from EXT unbiased file with emulated zero suppression
if extunbfile_emu:
    normfactor_emu_done = False
    for rootkey in extunbfile_emu.GetListOfKeys():
        if rootkey.GetClassName() != "TH1D":
            continue
        htemp = TH1D()
        htemp = rootkey.ReadObj()
        # Get normalization factor denominator
        if htemp.GetName().startswith("hitMultiplicity"):
            if normfactor_emu_done == False:
                normfactor_emu *= 1./(htemp.GetEntries()*extunbevtsize_emu)
                normfactor_emu_done = True
            # For hit multiplicity/event we have to correct for the different event size
            hnew = htemp.Clone()
            hnew.Reset()
            for ibin in xrange( htemp.GetNbinsX() ):
                # Trigger stream events are twice as long as SN stream ones, so they can have twice as many hits
                ibinnew = int( math.floor(ibin*snevtsize/extunbevtsize_emu) )
                hnew.SetBinContent( ibinnew, hnew.GetBinContent(ibinnew) + htemp.GetBinContent(ibin) )
            htemp = hnew    
        htemp.Sumw2()
        extunbhdict_emu.update( {htemp.GetName() : htemp} )

# Retrieve objects from EXT unbiased file with emulated zero suppression and hit-merger algorithm
if extunbfile_emumerged:
    normfactor_emumerged_done = False
    for rootkey in extunbfile_emumerged.GetListOfKeys():
        if rootkey.GetClassName() != "TH1D":
            continue
        htemp = TH1D()
        htemp = rootkey.ReadObj()
        # Get normalization factor denominator
        if htemp.GetName().startswith("Sum_hitMultiplicity"):
            if normfactor_emumerged_done == False:
                normfactor_emumerged *= 1./(htemp.GetEntries()*extunbevtsize_emu)
                normfactor_emumerged_done = True
            # For hit multiplicity/event we have to correct for the different event size
            hnew = htemp.Clone()
            hnew.Reset()
            for ibin in xrange( htemp.GetNbinsX() ):
                # Trigger stream events are twice as long as SN stream ones, so they can have twice as many hits
                ibinnew = int( math.floor(ibin*snevtsize/extunbevtsize_emu) )
                hnew.SetBinContent( ibinnew, hnew.GetBinContent(ibinnew) + htemp.GetBinContent(ibin) )
            htemp = hnew    
        htemp.Sumw2()
        extunbhdict_emumerged.update( {htemp.GetName() : htemp} )

# Turn off display of the EXT unbiased histograms
#normfactor[0] = 0.
#normfactor[1] = 0.
#normfactor[2] = 0.

# Lists to collect objects going out of scope
clist = []
extunbh_errorslist = []
leglist = []
hratiolist = []
line1list = []
extunbh_emu_collector = []
extunbh_emumerged_collector = []

# Loop for making plots
for snh in snhlist:

    c = TCanvas( "c_" + snh.GetName(), snh.GetTitle(), 800 if frontierpad else 1600, 800 )

    if frontierpad:
        ### Histograms plot ###
        pad1 = TPad( "pad1_" + snh.GetName(), snh.GetTitle(), 0., frontierpad, 1., 1. )
        pad1.SetBottomMargin(0) # Join with lower plot
        pad1.Draw()
        pad1.cd()

    extunbh = extunbhlist[snhlist.index(snh)]
    extunbh_integral = extunbh.Integral()
    # Apply normalization factor based on last character of name
    if extunbh.GetName().startswith("extunb_hEventHitMult"):
        # For the event hit multiplicity the event size is corrected differently
        extunbh.Scale( extunbevtsize/snevtsize*normfactor[int( snh.GetName()[-1] )] )
    else:
        extunbh.Scale( normfactor[int( snh.GetName()[-1] )] )
    extunbh.SetMarkerColor(2)
    extunbh.SetLineColor(2)
    extunbh_errors = extunbh.Clone( extunbh.GetName() + "_errors" )
    extunbh_errors.SetFillColor(46)
    extunbh_errors.Draw("E2") # Rectangle error boxes in the background
    #extunbh_errors.DrawNormalized("E2") # Rectangle error boxes in the background
    extunbh_errors.GetYaxis().SetRangeUser( 0.001, 
        1.1*extunbh_errors.GetMaximum() if ( extunbh_errors.GetMaximum() > snh.GetMaximum() ) else 1.1*snh.GetMaximum()  )
    #extunbh_errors.GetYaxis().SetLabelSize( (1. - frontierpad)*extunbh_errors.GetYaxis().GetLabelSize() )
    extunbh.Draw("hist same") # Central value line
    #extunbh.DrawNormalized("hist same") # Central value line

    #########################
    # Add emulated plots

    if len( extunbhdict_emu.keys() ):
        # Flipping-bits OFF
        extunbh_emu = None
        if extunbh.GetName().startswith("extunb_hEventHitMult"):
            if extunbh.GetName() == "extunb_hEventHitMult_0":
                extunbh_emu = extunbhdict_emu["hitMultiplicityfbOff_U"]
                print extunbh_emu.GetName()
            elif extunbh.GetName() == "extunb_hEventHitMult_1":
                extunbh_emu = extunbhdict_emu["hitMultiplicityfbOff_V"]
                print extunbh_emu.GetName()
            elif extunbh.GetName() == "extunb_hEventHitMult_2":
                extunbh_emu = extunbhdict_emu["hitMultiplicityfbOff_Y"]
                print extunbh_emu.GetName()
            extunbh_emu_integral = extunbh_emu.Integral()
            # For the event hit multiplicity the event size is corrected differently
            extunbh_emu.Scale( extunbevtsize_emu/snevtsize*normfactor_emu )
        elif extunbh.GetName().startswith("extunb_hEventHitSpectrum"):
            if extunbh.GetName() == "extunb_hEventHitSpectrum_0":
                extunbh_emu = extunbhdict_emu["hitIntegralfbOff_U"]
            elif extunbh.GetName() == "extunb_hEventHitSpectrum_1":
                extunbh_emu = extunbhdict_emu["hitIntegralfbOff_V"]
            elif extunbh.GetName() == "extunb_hEventHitSpectrum_2":
                extunbh_emu = extunbhdict_emu["hitIntegralfbOff_Y"]
            # Add rest of planes
            extunbh_emu_integral = extunbh_emu.Integral()
            extunbh_emu.Scale( normfactor_emu )
        if extunbh_emu:
            extunbh_emu.SetMarkerColor(3)
            extunbh_emu.SetLineColor(3)
            extunbh_emu_errors = extunbh_emu.Clone( extunbh_emu.GetName() + "_errors" )
            extunbh_emu_errors.SetFillColor(30)
            extunbh_emu_errors.SetFillStyle(3002)
            extunbh_emu_errors.Draw("E2 same") # Rectangle error boxes in the background
            #extunbh_emu_errors.DrawNormalized("E2 same") # Rectangle error boxes in the background
            extunbh_emu.Draw("hist same") # Central value line
            #extunbh.DrawNormalized("hist same") # Central value line
            # Keep the objects in memory
            extunbh_emu_collector.append( extunbh_emu )
            extunbh_emu_collector.append( extunbh_emu_errors )

            extunbh_errors.GetYaxis().SetRangeUser( 0.001, 
                                                    1.1*extunbh_errors.GetMaximum() if ( extunbh_errors.GetMaximum() > extunbh_emu.GetMaximum() ) else 1.1*extunbh_emu.GetMaximum()  )

        # Flipping-bits ON
        extunbh_emu = None
        if extunbh.GetName().startswith("extunb_hEventHitMult"):
            if extunbh.GetName() == "extunb_hEventHitMult_0":
                extunbh_emu = extunbhdict_emu["hitMultiplicityfbOn_U"]
                print extunbh_emu.GetName()
            elif extunbh.GetName() == "extunb_hEventHitMult_1":
                extunbh_emu = extunbhdict_emu["hitMultiplicityfbOn_V"]
                print extunbh_emu.GetName()
            elif extunbh.GetName() == "extunb_hEventHitMult_2":
                extunbh_emu = extunbhdict_emu["hitMultiplicityfbOn_Y"]
                print extunbh_emu.GetName()
            extunbh_emu_integral = extunbh_emu.Integral()
            # For the event hit multiplicity the event size is corrected differently
            extunbh_emu.Scale( extunbevtsize_emu/snevtsize*normfactor_emu )
        elif extunbh.GetName().startswith("extunb_hEventHitSpectrum"):
            if extunbh.GetName() == "extunb_hEventHitSpectrum_0":
                extunbh_emu = extunbhdict_emu["hitIntegralfbOn_U"]
            elif extunbh.GetName() == "extunb_hEventHitSpectrum_1":
                extunbh_emu = extunbhdict_emu["hitIntegralfbOn_V"]
            elif extunbh.GetName() == "extunb_hEventHitSpectrum_2":
                extunbh_emu = extunbhdict_emu["hitIntegralfbOn_Y"]
            # Add rest of planes
            extunbh_emu_integral = extunbh_emu.Integral()
            extunbh_emu.Scale( normfactor_emu )
        if extunbh_emu:
            extunbh_emu.SetMarkerColor(4)
            extunbh_emu.SetLineColor(4)
            extunbh_emu_errors = extunbh_emu.Clone( extunbh_emu.GetName() + "_errors" )
            extunbh_emu_errors.SetFillColor(40)
            extunbh_emu_errors.SetFillStyle(3003)
            extunbh_emu_errors.Draw("E2 same") # Rectangle error boxes in the background
            #extunbh_emu_errors.DrawNormalized("E2 same") # Rectangle error boxes in the background
            extunbh_emu.Draw("hist same") # Central value line
            #extunbh.DrawNormalized("hist same") # Central value line
            # Keep the objects in memory
            extunbh_emu_collector.append( extunbh_emu )
            extunbh_emu_collector.append( extunbh_emu_errors )

            extunbh_errors.GetYaxis().SetRangeUser( 0.001, 
                                                    1.1*extunbh_errors.GetMaximum() if ( extunbh_errors.GetMaximum() > extunbh_emu.GetMaximum() ) else 1.1*extunbh_emu.GetMaximum()  )

    if len( extunbhdict_emu.keys() ):
        # Flipping-bits ON + merged hits
        extunbh_emumerged = None
        if extunbh.GetName().startswith("extunb_hEventHitMult"):
            if extunbh.GetName() == "extunb_hEventHitMult_0":
                extunbh_emumerged = extunbhdict_emumerged["Sum_hitMultiplicityfbOn_U"]
                print extunbh_emumerged.GetName()
            elif extunbh.GetName() == "extunb_hEventHitMult_1":
                extunbh_emumerged = extunbhdict_emumerged["Sum_hitMultiplicityfbOn_V"]
                print extunbh_emumerged.GetName()
            elif extunbh.GetName() == "extunb_hEventHitMult_2":
                extunbh_emumerged = extunbhdict_emumerged["Sum_hitMultiplicityfbOn_Y"]
                print extunbh_emumerged.GetName()
            extunbh_emumerged_integral = extunbh_emumerged.Integral()
            # For the event hit multiplicity the event size is corrected differently
            extunbh_emumerged.Scale( extunbevtsize_emu/snevtsize*normfactor_emumerged )
            print normfactor_emumerged
        elif extunbh.GetName().startswith("extunb_hEventHitSpectrum"):
            if extunbh.GetName() == "extunb_hEventHitSpectrum_0":
                extunbh_emumerged = extunbhdict_emumerged["Sum_hitIntegralfbOn_U"]
            elif extunbh.GetName() == "extunb_hEventHitSpectrum_1":
                extunbh_emumerged = extunbhdict_emumerged["Sum_hitIntegralfbOn_V"]
            elif extunbh.GetName() == "extunb_hEventHitSpectrum_2":
                extunbh_emumerged = extunbhdict_emumerged["Sum_hitIntegralfbOn_Y"]
            # Add rest of planes
            extunbh_emumerged_integral = extunbh_emumerged.Integral()
            extunbh_emumerged.Scale( normfactor_emumerged )
            print normfactor_emumerged
        if extunbh_emumerged:
            extunbh_emumerged.SetMarkerColor(7)
            extunbh_emumerged.SetLineColor(7)
            extunbh_emumerged_errors = extunbh_emumerged.Clone( extunbh_emumerged.GetName() + "_errors" )
            extunbh_emumerged_errors.SetFillColor(33)
            extunbh_emumerged_errors.SetFillStyle(3144)
            extunbh_emumerged_errors.Draw("E2 same") # Rectangle error boxes in the background
            #extunbh_emumerged_errors.DrawNormalized("E2 same") # Rectangle error boxes in the background
            extunbh_emumerged.Draw("hist same") # Central value line
            #extunbh.DrawNormalized("hist same") # Central value line
            # Keep the objects in memory
            extunbh_emumerged_collector.append( extunbh_emumerged )
            extunbh_emumerged_collector.append( extunbh_emumerged_errors )

            extunbh_errors.GetYaxis().SetRangeUser( 0.001, 
                                                    1.1*extunbh_errors.GetMaximum() if ( extunbh_errors.GetMaximum() > extunbh_emumerged.GetMaximum() ) else 1.1*extunbh_emumerged.GetMaximum()  )

    # End of emulated plots
    #########################

    snh.SetMarkerStyle(21)
    snh.SetMarkerSize(0.5)
    snh.Draw("E P0 X0 same") # SN data
    #snh.DrawNormalized("E P0 X0 same") # SN data
    leg = TLegend(0.6, 0.7, 1.0, 1.0)
    leg.AddEntry( extunbh_errors, "Trigger stream (%i entries)" % int( extunbh_integral ), "FLP" )
    leg.AddEntry( snh, "SN stream (%i entries)" % int( snh.Integral() ), "LEP" )
    if logx: leg.SetFillStyle(0) # Make it transparent to mitigate the change of the distribution
    leg.Draw()
    if logx: gPad.SetLogx()
    if logy: 
        extunbh_errors.GetYaxis().SetRangeUser( 0.9, 1.1*extunbh_errors.GetMaximum() )
        gPad.SetLogy()
    gPad.Modified()
    gPad.Update()
    print "%s: rate %.2f +/- %.2f" % ( snh.GetName(), snh.Integral()/snexposure, math.sqrt(snh.Integral())/snexposure )
    print "\t trigger: rate %.2f +/- %.2f" % ( extunbh_integral/extunbexposure, math.sqrt(extunbh_integral)/extunbexposure )

    ### Ratio plot ###
    if frontierpad:
        c.cd()
        pad2 = TPad( "pad2_" + snh.GetName(), snh.GetTitle(), 0., 0.05, 1., frontierpad )
        pad2.SetTopMargin(0) # Join with upper plot
        pad2.SetBottomMargin(0.3)
        pad2.Draw()
        pad2.cd()

        hratio = snh.Clone()
        hratio.Sumw2()
        hratio.Divide(extunbh)
            
        hratio.GetYaxis().SetTitle("SN/Trigger")
        hratio.GetYaxis().SetRangeUser(0., 1.999)
        hratio.GetYaxis().SetNdivisions(505)
        hratio.GetYaxis().SetLabelSize( ((1. - frontierpad)/frontierpad)*extunbh_errors.GetYaxis().GetLabelSize() )
        hratio.GetYaxis().SetTitleSize( ((1. - frontierpad)/frontierpad)*extunbh_errors.GetYaxis().GetTitleSize() )
        hratio.GetYaxis().SetTitleOffset( frontierpad/(1. - frontierpad)*extunbh_errors.GetYaxis().GetTitleOffset() )
        hratio.GetYaxis().CenterTitle()
        hratio.GetXaxis().SetLabelSize( ((1. - frontierpad)/frontierpad)*extunbh_errors.GetXaxis().GetLabelSize() )
        hratio.GetXaxis().SetTitleSize( ((1. - frontierpad)/frontierpad)*extunbh_errors.GetXaxis().GetTitleSize() )
        hratio.SetMarkerStyle(21)
        hratio.SetMarkerSize(0.5)
        hratio.Draw("E P0 X0")
        line1 = TLine( extunbh_errors.GetXaxis().GetXmin(), 1., extunbh_errors.GetXaxis().GetXmax(), 1. )
        line1.SetLineWidth(1)
        line1.SetLineColor(2)
        line1.Draw()
        hratio.Draw("E P0 X0 same") # Draw again over line
        if logx: gPad.SetLogx()        
        # Keep the objects in memory
        line1list.append(line1)
        hratiolist.append(hratio)

    leglist.append(leg)
    extunbh_errorslist.append(extunbh_errors)
    c.Print( c.GetName() + ("_logx" if logx else "") + ("_logy" if logy else "") + ("_ratio.png" if frontierpad else ".png") )
    #c.Print( c.GetName() + ("_logx" if logx else "") + ("_bigratio.png" if frontierpad else ".png") )
    #c.Print( "norm_" + c.GetName() + ".png" )
    clist.append(c)

raw_input("Press the <ENTER> key to continue...")
