import ROOT
import numpy as np

avogadro_constant = 6.022140e23
# For a complete analysis there should be a cut on the vertex separation too, at the very least. These are
# sensitivity module branches

num_e = "sensitivity.number_of_electrons==2"
two_cal = "sensitivity.passes_two_calorimeters"
assoc_cal = "sensitivity.passes_associated_calorimeters"
ext_prob = "sensitivity.external_probability<0.01"
int_prob = "sensitivity.internal_probability>0.04"
two_clust = "sensitivity.passes_two_clusters"
two_track = "sensitivity.passes_two_tracks"
energy_cut = "(sensitivity.total_calorimeter_energy)>=2 && (sensitivity.total_calorimeter_energy)< 3.2"
main_cut_str = num_e + " && " + two_cal + " && " + assoc_cal + " && " + ext_prob + " && " + int_prob + " && " + two_clust + " && " + two_track

# check these numbers are the right way round
prob_cut_str = ""

selenium = {
    "name": "Se",
    "atomic_mass": 82,
    "2nbb_HL": 10.07e19,
    "min_energy": 2,  # MeV
    "max_energy": 3.2,  # MeV
    "mass": 7
}

thallium = {
    "name": "Tl",
    "atomic_mass": 208,
    "activity": 300 * 7  # µBq
}

thallium_foils_filename = "~/Desktop/sensitivity/thallium.root"

bismuth = {
    "name": "Bi",
    "atomic_mass": 214,
    "activity": 370  # µBq
}

bismuth_foils_filename = "~/Desktop/sensitivity/bismuth.root"


# This class is to set up a source isotope eg Se82. This allows you to calculate sensitivities for different detector
# configs
class Isotope:

    def __init__(self, description: dict):
        self.molar_mass_ = description["atomic_mass"]
        self.isotope_name_ = description["name"]

    def get_molar_mass(self):
        return self.molar_mass_

    def get_isotope_name(self):
        return self.isotope_name_

    def get_molar_mass_text(self):
        return str(self.molar_mass_)

    def set_isotope_name(self, val: str):
        self.isotope_name_ = val

    def set_molar_mass(self, val: int):
        self.molar_mass_ = val


# Background isotopes could be Bi214 or Tl208 in the foil or radon in the tracker
class BackgroundIsotope(Isotope):

    def __init__(self, description: dict, root_file_name: str, isotope_location: str):
        super().__init__(description)
        self.activity_ = description["activity"]
        self.root_file_name_ = root_file_name
        self.isotope_location_ = isotope_location

    def set_root_file_name(self, val: str):
        self.root_file_name_ = val

    def get_root_file_name(self):
        return self.root_file_name_

    def set_activity(self, val: float):
        self.activity_ = val

    def get_activity(self):
        return self.activity_

    def get_isotope_location(self):
        return self.isotope_location_

    def set_isotope_location(self, val: str):
        self.isotope_location_ = val


class SampleIsotope(Isotope):

    def __init__(self, description: dict,
                 exposure_years: float,
                 frac_events_2bbsample: float):
        super().__init__(description)
        self.bkgd_half_life_ = description["2nbb_HL"]
        self.min_energy_ = description["min_energy"]
        self.max_energy_ = description["max_energy"]
        self.isotope_mass_kg_ = description["mass"]
        self.frac_events_2bbsample_ = frac_events_2bbsample
        self.exposure_years_ = exposure_years

    def get_bkgd_half_life(self):
        return self.bkgd_half_life_

    def get_exposure_years(self):
        return self.exposure_years_

    def get_min_energy(self):
        return self.min_energy_

    def get_max_energy(self):
        return self.max_energy_

    def get_frac_events_2bbsample(self):
        return self.frac_events_2bbsample_

    def get_isotope_mass_kg(self):
        return self.isotope_mass_kg_

    def set_bkgd_half_life(self, val: float):
        self.bkgd_half_life_ = val

    def set_exposure_years(self, val: float):
        self.exposure_years_ = val

    def set_min_energy(self, val: float):
        self.min_energy_ = val

    def set_max_energy(self, val: float):
        self.max_energy_ = val

    def set_frac_events_2bbsample(self, val: float):
        self.frac_events_2bbsample_ = val

    def set_isotope_mass_kg(self, val: float):
        self.isotope_mass_kg_ = val

    def print_all(self):
        print("Isotope:", self.get_isotope_name(), "-", self.get_molar_mass())
        print("Mass:", self.get_isotope_mass_kg(), "kg")
        print("Exposure:", self.get_exposure_years(), "years")
        print("2nubb halflife:", self.get_bkgd_half_life(), "years")
        print("Energy range:", self.get_min_energy(), "-", self.get_max_energy(), "MeV contains",
              self.get_frac_events_2bbsample(), "of total events")


#  make a cut on true energy in the 2nubb simulation. That allows me to use a smaller simulated file
# that doesn't include events I would cut anyway. If you simulate like this, you can get the fraction of the total
# spectrum that was retained from the brio file. Put that number in here. If you simulate a full spectrum, it is 1.


def calculate_efficiencies(tree: ROOT.TTree, sample: SampleIsotope, is2nubb: bool, all_entries: int, identifier):
    print(">>> calculate_efficiencies()")
    is2nubb = False
    all_entries = 0

    type_ = sample.get_isotope_name() + sample.get_molar_mass_text() + "_" + identifier

    min_energy = 2
    max_energy = 3.6
    nbins = int(max_energy * 20 - min_energy * 20)

    # Reconstructable events
    tot_entries = tree.GetEntries()
    stack = ROOT.THStack("hs", "")

    # In the 2nubb case, some simulate events are not reconstructed
    print("Total events", tot_entries)
    print("Reconstructable entries", tot_entries)

    # Our efficiency should be calculated as a fraction of the total events we generated, not as a fraction of the events we managed to reconstruct

    if all_entries > 0:
        tot_entries = all_entries # Use this if not all entries were reconstructable

    passes_cal_cut = tree.GetEntries(two_cal)
    print("Two triggered calorimeters:", float(passes_cal_cut) / float(tot_entries) * 100, "%")
    name = "h_" + type_ + "_" + two_cal
    temp_0 = ROOT.TH1D(name, name, nbins, min_energy, max_energy)
    tree.Draw("(sensitivity.total_calorimeter_energy) >> " + name, two_cal, "HIST")
    temp_0.SetFillColor(7)
    temp_0.GetXaxis().SetTitle("E_{1}+E_{2} (MeV)")
    stack.Add(temp_0)
    '''c_0 = ROOT.TCanvas()
    temp_0.Draw()
    c_0.cd()
    c_0.SaveAs("~/Desktop/sensitivity/plots/" + type_ + "_0_0.png")'''


    # 2 tracker clusters
    passes_cluster_cut = tree.GetEntries(two_cal + " && " + two_clust)
    print("Two tracker clusters with a minimum of 3 cells:", float(passes_cluster_cut) / float(tot_entries) * 100, "%")
    name = "h_" + type_ + "_" + two_cal + "_" + two_clust
    temp_1 = ROOT.TH1D(name, name, nbins, min_energy, max_energy)
    tree.Draw("(sensitivity.total_calorimeter_energy) >> " + name, two_cal + " && " + two_clust, "HIST")
    temp_1.SetFillColor(1)
    stack.Add(temp_1)
    '''c_1 = ROOT.TCanvas()
    temp_1.Draw()
    c_1.cd()
    c_1.SaveAs("~/Desktop/sensitivity/plots/" + type_ + "_1_1.png")'''

    # 2 reconstructed tracks
    passes_tracker_cut = tree.GetEntries(two_cal + " && " + two_clust + " && " + two_track)
    print("Two reconstructed tracks:", float(passes_tracker_cut) / float(tot_entries) * 100, "%")
    name = "h_" + type_ + "_" + two_cal + "_" + two_clust + "_" + two_track
    temp_2 = ROOT.TH1D(name, name, nbins, min_energy, max_energy)
    tree.Draw("(sensitivity.total_calorimeter_energy) >> " + name, two_cal + " && " + two_clust + " && " + two_track, "HIST")
    temp_2.SetFillColor(2)
    stack.Add(temp_2)
    '''c_2 = ROOT.TCanvas()
    temp_2.Draw()
    c_2.cd()
    c_2.SaveAs("~/Desktop/sensitivity/plots/" + type_ + "_2_2.png")'''

    # Tracks have associated calorimeter hits
    passes_associated_cal_cut = tree.GetEntries(two_cal + " && " + two_clust + " && " + two_track + " && " + assoc_cal)
    print("Tracks have associated calorimeter hits:", float(passes_associated_cal_cut) / float(tot_entries) * 100, "%")
    name = "h_" + type_ + "_" + two_cal + "_" + two_clust + "_" + two_track + "_" + assoc_cal
    temp_3 = ROOT.TH1D(name, name, nbins, min_energy, max_energy)
    tree.Draw("(sensitivity.total_calorimeter_energy) >> " + name, two_cal + " && " + two_clust + " && " + two_track + " && " + assoc_cal, "HIST")
    temp_3.SetFillColor(3)
    stack.Add(temp_3)
    '''c_3 = ROOT.TCanvas()
    temp_3.Draw()
    c_3.cd()
    c_3.SaveAs("~/Desktop/sensitivity/plots/" + type_ + "_3_3.png")'''

    # Passes NEMO3 internal / external probability cuts
    passes_prob_cut = tree.GetEntries(two_cal + " && " + two_clust + " && " + two_track + " && " + assoc_cal + " && " +
                                      ext_prob + " && " + int_prob)
    print("And passes internal/external probability cut:", float(passes_prob_cut) / float(tot_entries) * 100, "%")
    name = "h_" + type_ + "_" + two_cal + "_" + two_clust + "_" + two_track + "_" + assoc_cal + "_ext_int_prob"
    temp_4 = ROOT.TH1D(name, name, nbins, min_energy, max_energy)
    tree.Draw("(sensitivity.total_calorimeter_energy) >> " + name, two_cal + " && " + two_clust + " && " + two_track + " && " + assoc_cal + " && " +
              ext_prob + " && " + int_prob, "HIST")
    temp_4.SetFillColor(4)
    stack.Add(temp_4)
    '''c_4 = ROOT.TCanvas()
    temp_4.Draw()
    c_4.cd()
    c_4.SaveAs("~/Desktop/sensitivity/plots/" + type_ + "_4_4.png")'''

    # Passes NEMO3 internal / external probability cuts, 2 - 3.2 MeV
    passes_prob_roi = tree.GetEntries(two_cal + " && " + two_clust + " && " + two_track + " && " + assoc_cal + " && " +
                                      ext_prob + " && " + int_prob + " && " + energy_cut)
    print("Passes internal/external probability in Se82ROI:", float(passes_prob_roi) / float(tot_entries) * 100, "%")
    name = "h_" + type_ + "_" + two_cal + "_" + two_clust + "_" + two_track + "_" + assoc_cal + "_ext_int_prob_energy_cut"
    temp_5 = ROOT.TH1D(name, name, nbins, min_energy, max_energy)
    tree.Draw("(sensitivity.total_calorimeter_energy) >> " + name, two_cal + " && " + two_clust + " && " + two_track + " && " + assoc_cal + " && " +
              ext_prob + " && " + int_prob + " && " + energy_cut, "HIST")
    temp_5.SetFillColor(5)
    stack.Add(temp_5)
    '''c_5 = ROOT.TCanvas()
    temp_5.Draw()
    c_5.cd()
    c_5.SaveAs("~/Desktop/sensitivity/plots/" + type_ + "_5_5.png")'''

    # And both tracks have negative charge reconstruction
    passes_two_electron_cut = tree.GetEntries(two_cal + " && " + two_clust + " && " + two_track + " && " + assoc_cal +
                                              " && " + ext_prob + " && " + int_prob + " && " + energy_cut + " && " + num_e)
    print("And both tracks have negative charge: ", float(passes_two_electron_cut) / float(tot_entries) * 100, "%")
    name = "h_" + type_ + "_" + two_cal + "_" + two_clust + "_" + two_track + "_" + assoc_cal + "_ext_int_prob_energy_cut_num_e"
    temp_6 = ROOT.TH1D(name, name, nbins, min_energy, max_energy)
    tree.Draw("(sensitivity.total_calorimeter_energy) >> " + name, two_cal + " && " + two_clust + " && " + two_track + " && " + assoc_cal + " && " +
              ext_prob + " && " + int_prob + " && " + energy_cut + " && " + num_e,  "HIST")
    temp_6.SetFillColor(6)
    stack.Add(temp_6)
    '''c_6 = ROOT.TCanvas()
    temp_6.Draw()
    c_6.cd()
    c_6.SaveAs("~/Desktop/sensitivity/plots/" + type_ + "_6_6.png")'''

    c = ROOT.TCanvas()
    c.cd()
    stack.Draw("nostack")
    stack.GetXaxis().SetTitle("E_{1}+E_{2} (MeV)")
    stack.GetYaxis().SetTitle("")
    stack.SetTitle("")
    c.SetLogy(True)
    stack_legend = ROOT.TLegend(0.2, 0.2, 0.35, 0.35)
    ROOT.gStyle.SetLegendBorderSize(0)
    stack_legend.AddEntry(temp_0, '2calos')
    stack_legend.AddEntry(temp_1, "2clusters")
    stack_legend.AddEntry(temp_2, '2tracks')
    stack_legend.AddEntry(temp_3, "associated calo")
    stack_legend.AddEntry(temp_4, 'ext and int probs')
    stack_legend.AddEntry(temp_5, "energy")
    stack_legend.AddEntry(temp_6, '2e-')
    stack_legend.Draw()
    c.SaveAs("~/Desktop/sensitivity/plots/" + type_ + "_stack.png")

    return tot_entries


def plot_bkgd_eff(bkgd_isotope: BackgroundIsotope, additional_cut: str, min_energy: float, max_energy: float):
    print(">>> plot_bkgd_eff()")
    file = ROOT.TFile(bkgd_isotope.get_root_file_name())
    #print(bkgd_isotope.get_root_file_name())
    tree = file.Get("Sensitivity;1")
    n_entries = tree.GetEntries()
    h_eff = plot_efficiency(tree, n_entries, additional_cut, min_energy, max_energy, "eff_"+bkgd_isotope.get_isotope_name())
    try:
        h_eff.GetEentries()
    except:
        print("None type plot_bkgd_eff")
    title = bkgd_isotope.get_isotope_name() + "-" + bkgd_isotope.get_molar_mass_text() + " (" + bkgd_isotope.get_isotope_location() + ")"
    h_eff.SetTitle(title)
    return h_eff


def plot_bkgd_energy(tree: ROOT.TTree, bkgd_isotope: BackgroundIsotope, additional_cut: str, min_energy: float, max_energy: float):
    print(">>> plot_bkgd_energy()")
    name = "energy" + bkgd_isotope.get_isotope_name()
    nbins = int(max_energy * 20 - min_energy * 20)
    h_energy = ROOT.TH1D(name, name, nbins, min_energy, max_energy)
    title = bkgd_isotope.get_isotope_name() + "-" + bkgd_isotope.get_molar_mass_text() + " (" + bkgd_isotope.get_isotope_location() + ")"
    h_energy.SetTitle(title)

    tree.Draw("(sensitivity.total_calorimeter_energy)>>"+name, main_cut_str + additional_cut, "HIST")
    '''tree.Draw("(sensitivity.total_calorimeter_energy)>>energy",
              "sensitivity.number_of_electrons==2 && sensitivity.passes_two_calorimeters && sensitivity.passes_associated_calorimeters " + additional_cut,
              "HIST")'''

    temp = ROOT.TCanvas()
    h_energy.GetXaxis().SetTitle("E_{1}+E_{2} (MeV)")
    h_energy.Draw()
    temp.SaveAs("~/Desktop/sensitivity/plots/" + name + ".png")
    del temp

    return h_energy


def plot_efficiency(tree: ROOT.TTree, total_entries: float, additional_cut: str, min_energy: float, max_energy: float, name: str):
    print(">>> plot_efficiency()")
    nbins = int((max_energy * 20) - (min_energy * 20))
    h_eff = ROOT.TH1D(name, name, nbins, min_energy, max_energy)

    for i_bin in range(nbins):
        low_energy_lim = h_eff.GetXaxis().GetBinLowEdge(i_bin)

        cut = main_cut_str + prob_cut_str + " &&(sensitivity.total_calorimeter_energy) >= {} && (sensitivity.total_calorimeter_energy) < {}".format(low_energy_lim, max_energy)
        cut += additional_cut

        entries_passing_cut = tree.GetEntries(cut)
        h_eff.SetBinContent(i_bin, float(entries_passing_cut) / float(total_entries))

    '''h_eff.SetMarkerSize(1)
    h_eff.SetMarkerStyle(10)'''
    h_eff.GetYaxis().SetTitle("Efficiency")
    title = 'Energy (MeV) $\leq$ $\Sigma_{12}$ $E_{calibrated}\leq$' + ' {:.1} MeV'.format(max_energy)
    h_eff.GetXaxis().SetTitle(title)

    temp = ROOT.TCanvas()
    h_eff.Draw()
    temp.SaveAs("~/Desktop/sensitivity/plots/"+name+".png")
    del temp

    return h_eff


def estimate_half_life_sensitivity(sig_eff: float, bkgd_eff: float, sample: SampleIsotope):
    print(">>> estimate_half_life_sensitivity()")
    num_of_sigma = 4.
    # Get the total number of background events we expect to see in the energy window, after all cuts
    bkgd_events = estimate_bkdg_events(bkgd_eff, sample)
    if bkgd_events == 0:
        return 0 # Is this true? Maybe if there is no background there is no background and that is just great

    # Use formula from p7 of Reviews of Modern Physics vol 80 issue 2 pages 481 - 516(2008)

    # Background per kilo per year = bDelta E in formula
    bkgd_rate = bkgd_events / (sample.get_isotope_mass_kg() * sample.get_exposure_years())

    sensitivity = (4.16e26 / num_of_sigma) * (sig_eff / sample.get_molar_mass()) * np.sqrt(sample.get_isotope_mass_kg() * sample.get_exposure_years() / bkgd_rate)

    return sensitivity


def estimate_bkdg_events(bkgd_eff: float, sample: SampleIsotope):
    print(">>> estimate_bkdg_events()")
    # Get the number of atoms you start with
    n_source_atoms = avogadro_constant * (sample.get_isotope_mass_kg() * 1000) / sample.get_molar_mass() # Avogadro is grams

    # Get the number of atoms there will be remaining after the exposure time in years, based on the halflife for 2nu double beta
    # n(t) = n(0) exp - (t / tau) where tau = halflife / ln2
    # double nRemainingAtoms = nSourceAtoms * TMath::Exp(-1 * TMath::Log(2) * exposureYears / backgroundHalfLife);

    # The exponent is around 10 ^ -20, it 's too small for TMath::Exp to deal with
    # Can we just go for a Taylor expansion for 1-e ^ -x where x is v small?
    # 1(- e ^ -x) ~ x so...

    tot_decays = n_source_atoms * (np.log(2) * sample.get_exposure_years() / sample.get_bkgd_half_life())

    # Multiply by the efficiency and that is the amount of background events you expect to see

    events = tot_decays * bkgd_eff

    print(tot_decays, "backgrounds, of which we see", events)
    return events


def isotope_plots(filename_0nbb: str, filename_2nbb: str, sample: SampleIsotope):
    print(">>> isotope_plots()")
    file_0nbb = ROOT.TFile(filename_0nbb)
    tree_0nbb = file_0nbb.Get("Sensitivity;1")

    file_2nbb = ROOT.TFile(filename_2nbb)
    tree_2nbb = file_2nbb.Get("Sensitivity;1")

    print("Sample:", sample.get_isotope_name(), "-", sample.get_molar_mass_text())
    print("---")
    print("0nubb sample:")
    print("---")

    tot_entries_0nbb = calculate_efficiencies(tree_0nbb, sample, False, 0, "0nbb")

    print("---")
    print("2nubb sample:")
    print("---")

    tot_entries_2nbb = calculate_efficiencies(tree_2nbb, sample, True, 0, "2nbb")

    ROOT.gStyle.SetOptStat(0)

    '''extra_cut_plots(tree_0nbb, tot_entries_0nbb, tree_2nbb, tot_entries_2nbb,
                        "&& sensitivity.number_of_electrons==2", "Two electron tracks required", "two_electron_cut",
                        sample)'''
    '''extra_cut_plots(tree_0nbb, tot_entries_0nbb, tree_2nbb, tot_entries_2nbb,
                        "&& sensitivity.number_of_electrons>=1", "One or more negatively charged tracks required",
                        "one_electron_cut", sample)'''

    extra_cut_plots(tree_0nbb, tot_entries_0nbb, tree_2nbb, tot_entries_2nbb, "",
                    "", "", sample)

    '''extra_cut_plots(tree_0nbb, tot_entries_0nbb, tree_2nbb, tot_entries_2nbb, "", "No probability cuts",
                    "no_electron_cut", sample)'''

    '''extra_cut_plots(tree_0nbb, tot_entries_0nbb, tree_2nbb, tot_entries_2nbb,
                        "&& sensitivity.number_of_electrons==2 && sensitivity.external_probability<0.01 && sensitivity.internal_probability>0.04",
                        "Two -ve tracks,  NEMO3 internal/external prob cuts", "two_electron_int_ext_prob", sample)'''


def extra_cut_plots(tree_0nbb: ROOT.TTree, tot_entries_0nbb: float, tree_2nbb: ROOT.TTree, tot_entries_2nbb: float,
                    extra_cut: str, cut_title: str, cut_filename_suffix: str, sample: SampleIsotope):
    print(">>> extra_cut_plots()")
    nbins = int((sample.get_max_energy() * 20 - sample.get_min_energy() * 20))

    c = ROOT.TCanvas()
    c.cd()

    h_energy_0nbb = ROOT.TH1D("h_energy_0nbb", "h_energy_0nbb", nbins, sample.get_min_energy(), sample.get_max_energy())
    h_energy_0nbb.SetLineColor(1)
    h_energy_0nbb.GetXaxis().SetTitle("E_{1}+E_{2} (MeV)")
    h_energy_0nbb.GetYaxis().SetTitle('Counts - normalised $0\\nu\\beta\\beta$')
    h_energy_0nbb.SetTitle("^{" + sample.get_molar_mass_text() + "}" + sample.get_isotope_name() + " " + cut_title)

    total_cut = main_cut_str + extra_cut
    print("Cut is: ", total_cut)

    tree_0nbb.Draw("(sensitivity.total_calorimeter_energy) >> h_energy_0nbb", total_cut, "HIST")
    '''tree_0nbb.Draw("(sensitivity.total_calorimeter_energy)>>energy0nuBB",
                   "sensitivity.number_of_electrons==2 && sensitivity.passes_two_calorimeters && sensitivity.passes_associated_calorimeters && (sensitivity.higher_electron_energy != sensitivity.lower_electron_energy) "
                    + extra_cut,
                   "HIST")'''

    h_energy_2nbb = ROOT.TH1D("h_energy_2nbb", "h_energy_2nbb", nbins, sample.get_min_energy(), sample.get_max_energy())
    h_energy_2nbb.SetLineColor(2)

    tree_2nbb.Draw("(sensitivity.total_calorimeter_energy) >> h_energy_2nbb", total_cut, "HIST SAME")

    h_energy_2nbb.Scale(h_energy_0nbb.Integral() / h_energy_2nbb.Integral())

    energy_legend = ROOT.TLegend(0.2, 0.2, 0.35, 0.35)
    ROOT.gStyle.SetLegendBorderSize(0)
    energy_legend.AddEntry(h_energy_0nbb, '$0\\nu\\beta\\beta$')
    energy_legend.AddEntry(h_energy_2nbb, "$2\\nu\\beta\\beta$")
    c.SetLogy(True)
    energy_legend.Draw()

    c.SaveAs("~/Desktop/sensitivity/plots/" + sample.get_isotope_name() +
             sample.get_molar_mass_text() + "_energy_" + cut_filename_suffix + ".png")
    del c

    h_eff_0nbb = plot_efficiency(tree_0nbb, tot_entries_0nbb, extra_cut, sample.get_min_energy(),
                                 sample.get_max_energy(), "eff_0nbb")
    h_eff_2nbb = plot_efficiency(tree_2nbb, tot_entries_2nbb, extra_cut, sample.get_min_energy(),
                                 sample.get_max_energy(), "eff_2nbb")

    # Plot efficiency for 0nbb and 2nbb
    h_eff_0nbb.SetLineColor(1)
    h_eff_2nbb.SetLineColor(2)

    h_eff_0nbb.SetTitle("^{" + sample.get_molar_mass_text() + "}" + sample.get_isotope_name() + " " + cut_title)
    #h_eff_0nbb.GetYaxis().SetRangeUser(1e-8, 1)

    c = ROOT.TCanvas()
    c.cd()

    h_eff_0nbb.Draw("HIST")
    h_eff_2nbb.Scale(sample.get_frac_events_2bbsample())
    h_eff_2nbb.Draw("HIST SAME")

    eff_legend = ROOT.TLegend(0.2, 0.2, 0.35, 0.35)
    ROOT.gStyle.SetLegendBorderSize(0)
    eff_legend.AddEntry(h_eff_0nbb, '$0\\nu\\beta\\beta$')
    eff_legend.AddEntry(h_eff_2nbb, "$2\\nu\\beta\\beta$")
    c.SetLogy(True)
    eff_legend.Draw()
    c.SaveAs("~/Desktop/sensitivity/plots/" + sample.get_isotope_name() +
             sample.get_molar_mass_text() + "_efficiency_" + cut_filename_suffix + ".png")
    del c

    # Calculate efficiencies for other background isotopes
    bkgd_isotope_effs = ROOT.TList()
    bkgd_isotope_energies = ROOT.TList()
    bkgd_isotopes = []

    # Background isotopes go here
    # Newest measurement 370 uBq total

    tl_bkgd = BackgroundIsotope(thallium, thallium_foils_filename, "foils")
    bi_bkgd = BackgroundIsotope(bismuth, bismuth_foils_filename, "foils")

    bkgd_isotopes.append(tl_bkgd)
    bkgd_isotopes.append(bi_bkgd)

    '''# plot the bkgd eff Tl
    tl_file = ROOT.TFile(tl_bkgd.get_root_file_name())
    tl_tree = tl_file.Get("Sensitivity;1")
    n_entries = tl_tree.GetEntries()
    name = "eff_" + tl_bkgd.get_isotope_name()

    n_bins = int((sample.get_max_energy() * 20) - (sample.get_min_energy() * 20))
    h_tl_eff = ROOT.TH1D(name, name, nbins, sample.get_min_energy(), sample.get_max_energy())

    try:
        h_tl_eff.GetEntries()
    except:
        print("EFF 1 Null objects created")
        exit()

    for i_bin in range(nbins):
        low_energy_lim = h_tl_eff.GetXaxis().GetBinLowEdge(i_bin)

        cut = main_cut_str + prob_cut_str + \
              " &&(sensitivity.total_calorimeter_energy) >= {} && (sensitivity.total_calorimeter_energy) < {}".format(
                  low_energy_lim, sample.get_max_energy())
        cut += extra_cut

        entries_passing_cut = tl_tree.GetEntries(cut)
        h_tl_eff.SetBinContent(i_bin, float(entries_passing_cut) / float(n_entries))

    h_tl_eff.SetMarkerSize(1)
    h_tl_eff.SetMarkerStyle(10)
    h_tl_eff.GetYaxis().SetTitle("Efficiency")
    x_title = 'Energy (MeV) $\leq$ $\Sigma_{12}$ $E_{calibrated}\leq$' + ' {:.1} MeV'.format(sample.get_max_energy())
    h_tl_eff.GetXaxis().SetTitle(x_title)
    temp_canvas = ROOT.TCanvas()
    h_tl_eff.Draw()
    temp_canvas.SaveAs("~/Desktop/sensitivity/plots/" + name + ".png")
    del temp_canvas
    title = tl_bkgd.get_isotope_name() + "-" + tl_bkgd.get_molar_mass_text() + " (" + tl_bkgd.get_isotope_location() + ")"
    h_tl_eff.SetTitle(title)

    try:
        h_tl_eff.GetEntries()
    except:
        print("EFF 2 Null objects created")
        exit()

    # plot the bkgd energy
    temp_name = "energy" + tl_bkgd.get_isotope_name()
    h_tl_energy = ROOT.TH1D(temp_name, temp_name, n_bins, sample.get_min_energy(), sample.get_max_energy())
    try:
        h_tl_energy.GetEntries()
    except:
        print("Energy 1 Null objects created")
        exit()

    title = tl_bkgd.get_isotope_name() + "-" + tl_bkgd.get_molar_mass_text() + " (" + tl_bkgd.get_isotope_location() + ")"
    h_tl_energy.SetTitle(title)
    tl_tree.Draw("(sensitivity.total_calorimeter_energy)>>" + temp_name, main_cut_str + extra_cut, "HIST")
    temp_canvas = ROOT.TCanvas()
    h_tl_energy.Draw()
    temp_canvas.SaveAs("~/Desktop/sensitivity/plots/" + temp_name + ".png")
    del temp_canvas

    try:
        h_tl_energy.GetEntries()
    except:
        print("Energy 2 Null objects created")
        exit()

    # plot the bkgd eff Bi
    bi_file = ROOT.TFile(bi_bkgd.get_root_file_name())
    bi_tree = bi_file.Get("Sensitivity;1")
    n_entries = bi_tree.GetEntries()
    name = "eff_" + bi_bkgd.get_isotope_name()

    n_bins = int((sample.get_max_energy() * 20) - (sample.get_min_energy() * 20))
    h_bi_eff = ROOT.TH1D(name, name, nbins, sample.get_min_energy(), sample.get_max_energy())

    try:
        h_bi_eff.GetEntries()
    except:
        print("EFF 1 Null objects created")
        exit()

    for i_bin in range(nbins):
        low_energy_lim = h_bi_eff.GetXaxis().GetBinLowEdge(i_bin)

        cut = main_cut_str + prob_cut_str + \
              " &&(sensitivity.total_calorimeter_energy) >= {} && (sensitivity.total_calorimeter_energy) < {}".format(
                  low_energy_lim, sample.get_max_energy())
        cut += extra_cut

        entries_passing_cut = bi_tree.GetEntries(cut)
        h_bi_eff.SetBinContent(i_bin, float(entries_passing_cut) / float(n_entries))

    h_bi_eff.SetMarkerSize(1)
    h_bi_eff.SetMarkerStyle(10)
    h_bi_eff.GetYaxis().SetTitle("Efficiency")
    x_title = 'Energy (MeV) $\leq$ $\Sigma_{12}$ $E_{calibrated}\leq$' + ' {:.1} MeV'.format(sample.get_max_energy())
    h_bi_eff.GetXaxis().SetTitle(x_title)
    temp_canvas = ROOT.TCanvas()
    h_bi_eff.Draw()
    temp_canvas.SaveAs("~/Desktop/sensitivity/plots/" + name + ".png")
    del temp_canvas
    title = bi_bkgd.get_isotope_name() + "-" + bi_bkgd.get_molar_mass_text() + " (" + bi_bkgd.get_isotope_location() + ")"
    h_bi_eff.SetTitle(title)

    try:
        h_bi_eff.GetEntries()
    except:
        print("EFF 2 Null objects created")
        exit()

    # plot the bkgd energy
    temp_name = "energy" + tl_bkgd.get_isotope_name()
    h_bi_energy = ROOT.TH1D(temp_name, temp_name, n_bins, sample.get_min_energy(), sample.get_max_energy())
    try:
        h_bi_energy.GetEntries()
    except:
        print("Energy 1 Null objects created")
        exit()

    t =8
    title = bi_bkgd.get_isotope_name() + "-" + bi_bkgd.get_molar_mass_text() + " (" + bi_bkgd.get_isotope_location() + ")"
    h_bi_energy.SetTitle(title)
    bi_tree.Draw("(sensitivity.total_calorimeter_energy)>>" + temp_name, main_cut_str + extra_cut, "HIST")
    temp_canvas = ROOT.TCanvas()
    h_bi_energy.Draw()
    temp_canvas.SaveAs("~/Desktop/sensitivity/plots/" + temp_name + ".png")
    del temp_canvas

    try:
        h_bi_energy.GetEntries()
    except:
        print("Energy 2 Null objects created")
        exit()'''

    '''for i in range(len(bkgd_isotopes)):
        print("i:",i)
        # plot the bkgd eff
        temp_file = ROOT.TFile(bkgd_isotopes[i].get_root_file_name())
        temp_tree = temp_file.Get("Sensitivity;1")
        n_entries = temp_tree.GetEntries()
        name = "eff_" + bkgd_isotopes[i].get_isotope_name()

        n_bins = int((sample.get_max_energy() * 20) - (sample.get_min_energy() * 20))
        h_temp_eff = ROOT.TH1D(name, name, nbins, sample.get_min_energy(), sample.get_max_energy())

        try:
            h_temp_eff.GetEntries()
        except:
            print("EFF 1 Null objects created")
            exit()

        for i_bin in range(nbins):
            low_energy_lim = h_temp_eff.GetXaxis().GetBinLowEdge(i_bin)

            cut = main_cut_str + prob_cut_str + \
                  " &&(sensitivity.total_calorimeter_energy) >= {} && (sensitivity.total_calorimeter_energy) < {}".format(
                      low_energy_lim, sample.get_max_energy())
            cut += extra_cut

            entries_passing_cut = temp_tree.GetEntries(cut)
            h_temp_eff.SetBinContent(i_bin, float(entries_passing_cut) / float(n_entries))

        h_temp_eff.SetMarkerSize(1)
        h_temp_eff.SetMarkerStyle(10)
        h_temp_eff.GetYaxis().SetTitle("Efficiency")
        x_title = 'Energy (MeV) $\leq$ $\Sigma_{12}$ $E_{calibrated}\leq$' + ' {:.1} MeV'.format(sample.get_max_energy())
        h_temp_eff.GetXaxis().SetTitle(x_title)
        temp_canvas = ROOT.TCanvas()
        h_temp_eff.Draw()
        temp_canvas.SaveAs("~/Desktop/sensitivity/plots/" + name + ".png")
        del temp_canvas
        title = bkgd_isotopes[i].get_isotope_name() + "-" + bkgd_isotopes[i].get_molar_mass_text() + " (" + bkgd_isotopes[i].get_isotope_location() + ")"
        h_temp_eff.SetTitle(title)

        try:
            h_temp_eff.GetEntries()
        except:
            print("EFF 2 Null objects created")
            exit()

        # plot the bkgd energy
        temp_name = "energy" + bkgd_isotopes[i].get_isotope_name()
        h_temp_energy = ROOT.TH1D(temp_name, temp_name, n_bins, sample.get_min_energy(), sample.get_max_energy())
        try:
            h_temp_energy.GetEntries()
        except:
            print("Energy 1 Null objects created")
            exit()

        title = bkgd_isotopes[i].get_isotope_name() + "-" + bkgd_isotopes[i].get_molar_mass_text() + " (" + bkgd_isotopes[i].get_isotope_location() + ")"
        h_temp_energy.SetTitle(title)
        temp_tree.Draw("(sensitivity.total_calorimeter_energy)>>" + temp_name, main_cut_str + extra_cut, "HIST")
        temp_canvas = ROOT.TCanvas()
        h_temp_energy.Draw()
        temp_canvas.SaveAs("~/Desktop/sensitivity/plots/" + temp_name + ".png")
        del temp_canvas

        try:
            h_temp_energy.GetEntries()
        except:
            print("Energy 2 Null objects created")
            exit()

        try:
            h_temp_eff.GetEntries()
            h_temp_energy.GetEntries()
        except:
            print("One of Null objects created")
            exit()

        bkgd_isotope_effs.append(h_temp_eff)
        bkgd_isotope_energies.append(h_temp_energy)'''

    for i in range(len(bkgd_isotopes)):
        print(">>> Isotope:", bkgd_isotopes[i].get_isotope_name())
        file = ROOT.TFile(bkgd_isotopes[i].get_root_file_name())
        tree = file.Get("Sensitivity;1")
        n_entries = tree.GetEntries()

        bkgd_h_eff = plot_efficiency(tree, n_entries, extra_cut, sample.get_min_energy(), sample.get_max_energy(),
                                     "eff_" + bkgd_isotopes[i].get_isotope_name())
        title = bkgd_isotopes[i].get_isotope_name() + "-" + bkgd_isotopes[i].get_molar_mass_text() + " (" + bkgd_isotopes[i].get_isotope_location() + ")"
        bkgd_h_eff.SetTitle(title)

        bkgd_h_energy = plot_bkgd_energy(tree, bkgd_isotopes[i], extra_cut, sample.get_min_energy(), sample.get_max_energy())

        bkgd_isotope_effs.Add(bkgd_h_eff)
        bkgd_isotope_energies.Add(bkgd_h_energy)

    # Use the efficiencies to calculate a sensitivity
    c = ROOT.TCanvas()
    c.cd()
    h_sensitivity = estimate_sensitivity(h_energy_0nbb, h_energy_2nbb, sample, h_eff_0nbb, h_eff_2nbb)
    c.SetLogy(False)
    h_sensitivity.GetYaxis().SetRangeUser(0, 1e25)
    h_sensitivity.SetTitle("^{" + sample.get_molar_mass_text() + "}" + sample.get_isotope_name() + " " + cut_title)
    h_sensitivity.Draw("HIST")
    c.SaveAs("~/Desktop/sensitivity/plots/" + sample.get_isotope_name() +
             sample.get_molar_mass_text() + "_sensitivity_window_method_" + cut_filename_suffix + ".png")
    del c

    c = ROOT.TCanvas()
    c.cd()

    # Get an overall sensitivity from TLimit and compare
    h_temp_signal = h_energy_0nbb.Clone()
    h_scaled_bkgd = h_energy_2nbb.Clone()

    # Expected number of background events in the plotted range double
    bkgd_events = estimate_bkdg_events(h_eff_2nbb.GetBinContent(1), sample)

    print("Expected number of 2nubb events:", bkgd_events)
    h_scaled_bkgd.Sumw2()

    # If we scale by this number for our initial plot, we have it correctly-normalized
    h_scaled_bkgd.Scale(bkgd_events / h_scaled_bkgd.Integral())
    h_temp_data = h_scaled_bkgd.Clone()

    exp_sig_event_lim = exp_limit_sig_events(0.1, h_temp_signal, h_scaled_bkgd, h_temp_data)
    tlim_sensitivity = h_eff_0nbb.GetBinContent(1) * sample.get_isotope_mass_kg() * 1000 * \
                       sample.get_exposure_years() * avogadro_constant * \
                       np.log(2) / (sample.get_molar_mass() * exp_sig_event_lim)

    # Plot them all together
    # energy 0nubb scaled to...expected number of events or something?
    # energy of 2nubb scaled
    # energy of each background scaled

    h_scaled_bkgd.SetLineColor(4)
    h_scaled_bkgd.SetLineWidth(3)
    h_scaled_bkgd.GetYaxis().SetTitle("Expected events")
    h_scaled_bkgd.GetXaxis().SetTitle("Summed electron energies (MeV)")
    h_scaled_bkgd.SetTitle("Signal and background energy spectra")
    h_scaled_bkgd.Draw("HIST")
    legend = ROOT.TLegend(0.6, 0.65, 0.9, 0.85)
    ROOT.gStyle.SetLegendBorderSize(0)
    legend.AddEntry(h_scaled_bkgd, '$2\\nu\\beta\\beta$', "l")

    # Get TLimit sensitivity with additional backgrounds
    h_tot_bkgd = h_scaled_bkgd.Clone()
    h_tot_bkgd.Sumw2()
    colours = [7, 5, 6, 8, 9]

    for i in range(len(bkgd_isotopes)):
        this_bkgd_iso = bkgd_isotopes[i]
        h_this_iso_eff = bkgd_isotope_effs[i]
        h_this_iso_energy = bkgd_isotope_energies[i]
        '''if i == 0:
            this_bkgd_iso = tl_bkgd
            h_this_iso_eff = h_tl_eff.GetBinContent(1)
            h_this_iso_energy = h_tl_energy
        else: # i == 1:
            this_bkgd_iso = bi_bkgd
            h_this_iso_eff = h_bi_eff.GetBinContent(1)
            h_this_iso_energy = h_bi_energy'''

        this_iso_events = this_bkgd_iso.get_activity() * 1e-6 * h_this_iso_eff * (sample.get_exposure_years() * 3600 * 24 * 365.25)
        # Number of events; get the units right

        print(this_bkgd_iso.get_isotope_name(), this_bkgd_iso.get_isotope_location(), "Efficiency:", h_this_iso_eff,
              "Events:", this_iso_events)

        h_this_iso_energy.Sumw2()
        h_this_iso_energy.Scale(this_iso_events / h_this_iso_energy.Integral())

        h_tot_bkgd.Add(h_this_iso_energy)
        h_this_iso_energy.SetLineColor(colours[i])
        h_this_iso_energy.SetLineWidth(3)
        h_this_iso_energy.Draw("HIST SAME")
        legend.AddEntry(h_this_iso_energy,
                        "^{" + this_bkgd_iso.get_molar_mass_text() + "}" + this_bkgd_iso.get_isotope_name() + " (" + this_bkgd_iso.get_isotope_location() + ")",
                        "l")

    h_temp_data = h_tot_bkgd.Clone()
    tot_exp_sig_event_lim = exp_limit_sig_events(0.1, h_temp_signal, h_tot_bkgd, h_temp_data)
    tot_tlim_sensitivity = h_eff_0nbb.GetBinContent(
        1) * sample.get_isotope_mass_kg() * 1000 * sample.get_exposure_years() * avogadro_constant * np.log(2) / (
                                       sample.get_molar_mass() * tot_exp_sig_event_lim)

    # Print the results last because TLimit has masse of annoying output
    print("Sensitivity from TLimit (no other isotopes):", tlim_sensitivity, " years ")
    '''print("Sensitivity from Window Method (no other isotopes):",
          h_sensitivity.GetBinContent(h_sensitivity.GetMaximumBin()), "cutting at",
          h_sensitivity.GetBinLowEdge(h_sensitivity.GetMaximumBin()), "MeV")'''
    print("Sensitivity from TLimit including background isotopes:", tot_tlim_sensitivity, "years")

    h_scaled_sig = h_energy_0nbb.Clone()
    h_scaled_sig.SetLineColor(2)
    h_scaled_sig.SetLineWidth(3)
    print("Expected signal event limit:", tot_exp_sig_event_lim)
    h_scaled_sig.Scale(tot_exp_sig_event_lim/h_energy_0nbb.Integral())
    h_scaled_sig.Draw("HIST SAME")
    legend.AddEntry(h_scaled_sig, '$0\\nu\\beta\\beta$', "l")
    legend.Draw()
    c.SetLogy(True)
    c.SaveAs("~/Desktop/sensitivity/plots/all_bkgds.png")
    del c


def exp_limit_sig_events(conf_level: float, h_signal: ROOT.TH1D, h_bkgd: ROOT.TH1D, h_data: ROOT.TH1D):
    print(">>> exp_limit_sig_events()")
    low_bound = 0.1 / h_signal.Integral()
    high_bound = 1000.0 / h_signal.Integral()

    null_hyp_signal = h_signal.Clone("null_hyp_signal")
    null_hyp_signal.Scale(low_bound)
    disc_hyp_signal = h_signal.Clone("disc_hyp_signal")
    disc_hyp_signal.Scale(high_bound)

    data_source = ROOT.TLimitDataSource(null_hyp_signal, h_bkgd, h_data)

    confidence = ROOT.TLimit.ComputeLimit(data_source, 50000)
    low_bound_cl = confidence.CLs()

    del data_source

    data_source = ROOT.TLimitDataSource(disc_hyp_signal, h_bkgd, h_data)

    confidence = ROOT.TLimit.ComputeLimit(data_source, 50000)

    high_bound_cl = confidence.CLs()

    del data_source

    accuracy = 0.01
    this_cl = 0
    this_val = 0

    while abs(high_bound - low_bound) * h_signal.Integral() > accuracy:
        # bisection
        this_val = low_bound+(high_bound - low_bound) / 3

        h_this_signal = h_signal.Clone("test_signal")
        h_this_signal.Scale(this_val)

        data_source = ROOT.TLimitDataSource(h_this_signal, h_bkgd, h_data)
        confidence = ROOT.TLimit.ComputeLimit(data_source, 50000)

        this_cl = confidence.GetExpectedCLs_b()
        if this_cl > conf_level:
            low_bound = this_val
            low_bound_cl = this_cl
        else:
            high_bound = this_val
            high_bound_cl = this_cl


        del data_source
        del h_this_signal
        del confidence


    del null_hyp_signal
    del disc_hyp_signal

    return h_signal.Integral() * this_val


def estimate_sensitivity(h_energy_0nbb: ROOT.TH1D, h_energy_2nbb: ROOT.TH1D, sample: SampleIsotope,
                         eff_0nbb: ROOT.TH1D, eff_2nbb: ROOT.TH1D):
    print(">>> estimate_sensitivity()")
    h_scaled_bkgd = h_energy_2nbb.Clone()
    bkgd_events = estimate_bkdg_events(eff_2nbb.GetBinContent(1), sample)
    h_scaled_bkgd.Sumw2()
    # If we scale by this number for our initial plot, we have a correctly-normalized
    h_scaled_bkgd.Scale(bkgd_events / h_scaled_bkgd.Integral())
    h_sensitivity = h_energy_0nbb.Clone()
    h_sensitivity.Reset()

    h_sig_event_lim = h_energy_0nbb.Clone()
    h_sig_event_lim.Reset()

    print("Window method")
    for i in range(h_energy_0nbb.GetNbinsX()):
        if int(h_scaled_bkgd.Integral(i, h_energy_0nbb.GetNbinsX())) > 0:
            exp_sig_event_lim = window_method_find_exp_sig_evts(h_scaled_bkgd.Integral(i, h_energy_0nbb.GetNbinsX()))
            print("Bin", i, ": expected signal event limit:", exp_sig_event_lim)
            this_sensitivity = eff_0nbb.GetBinContent(i) * sample.get_isotope_mass_kg() * 1000 *sample.get_exposure_years() * avogadro_constant * np.log(2) / (sample.get_molar_mass() * exp_sig_event_lim)
            h_sensitivity.SetBinContent(i, this_sensitivity)
            h_sig_event_lim.SetBinContent(i, exp_sig_event_lim)

    h_sensitivity.SetMarkerSize(1)
    h_sensitivity.SetMarkerStyle(10)
    h_sensitivity.GetYaxis().SetTitle('$0\\nu\\beta\\beta$ halflife sensitivity')
    h_sensitivity.GetYaxis().SetTitleOffset(1.2)
    title = 'Energy (MeV) $\leq \Sigma_{12} E_{calibrated}\leq' + '{:.1f} MeV'.format(sample.get_max_energy())
    h_sensitivity.GetXaxis().SetTitle(title)

    h_sig_event_lim.Print("ALL")
    return h_sensitivity


def window_method_find_exp_sig_evts(B: float):
    print(">>> window_method_find_exp_sig_evts()")
    # Find S using CL(S + B) / CL(B) = 0.1 with N_Obs = B
    likelihood = 1.
    S = 0
    n_events = int(B)

    while likelihood > 0.999:
        S += 0.001
        CLsb = 0
        CLb = 0
        for i in range(n_events):
            CLsb += ROOT.TMath.Poisson(i, S + B)
            CLb += ROOT.TMath.Poisson(i, B)

        likelihood = CLsb / CLb
        #print("S =", S, "tCLb =", CLb, "tCLsb =", CLsb, "tCLs =", likelihood)

    return S


def main():
    print(">>> main()")
    exposure = 2.5  # years
    fraction = 0.04  # fraction of the full spectrum in our 2MeV+ sample

    se_sample = SampleIsotope(selenium, exposure, fraction)
    isotope_plots("~/Desktop/sensitivity/selenium_0nbb.root", "~/Desktop/sensitivity/selenium_2nbb.root", se_sample)


if __name__ == '__main__':
    main()
