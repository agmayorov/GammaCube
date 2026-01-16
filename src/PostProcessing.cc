#include "PostProcessing.hh"

namespace fs = std::filesystem;

PostProcessing::PostProcessing(std::string rootFilePath,
                               std::string outputFolderName,
                               bool processSecondaries,
                               bool processOptics,
                               double eMinMeV,
                               double eMaxMeV,
                               double eTrigMeV,
                               std::string particle)
    : rootFilePath(std::move(rootFilePath)),
      outputFolderName(std::move(outputFolderName)),
      processSecondaries(processSecondaries),
      processOptics(processOptics),
      eMinMeV(eMinMeV),
      eMaxMeV(eMaxMeV),
      eTrigMeV(eTrigMeV),
      particleName(std::move(particle)) {
    ROOT::EnableImplicitMT();
    gROOT->SetBatch(kTRUE);

    gErrorIgnoreLevel = kWarning;

    OpenRootFile();
    PrepareOutputDirs();
}

PostProcessing::~PostProcessing() = default;

void PostProcessing::OpenRootFile() {
    rootFile.reset(TFile::Open(rootFilePath.c_str(), "READ"));
    if (!rootFile || rootFile->IsZombie()) {
        throw std::runtime_error("Failed to open ROOT file: " + rootFilePath);
    }
}

void PostProcessing::PrepareOutputDirs() {
    fs::path rootPath(rootFilePath);
    fs::path rootDir = rootPath.parent_path();

    postProcessingDir = (rootDir / "post_processing").string();
    runDir = (fs::path(postProcessingDir) / outputFolderName).string();

    effectiveAreaDir = (fs::path(runDir) / "effective_area").string();
    sensitivityDir = (fs::path(runDir) / "sensitivity").string();
    csvDir = (fs::path(runDir) / "CSV").string();
    histogramsDir = (fs::path(runDir) / "Histograms").string();


    fs::create_directories(postProcessingDir);


    if (fs::exists(runDir)) {
        fs::remove_all(runDir);
    }

    fs::create_directories(effectiveAreaDir);
    fs::create_directories(sensitivityDir);
    fs::create_directories(csvDir);
    fs::create_directories(histogramsDir);
}

TH1* PostProcessing::GetHistOrThrow(const std::string& histName) {
    TH1* h = nullptr;
    rootFile->GetObject(histName.c_str(), h);
    if (!h) {
        throw std::runtime_error("Histogram not found: " + histName);
    }
    return h;
}

double PostProcessing::GeomCenter(double eLow, double eHigh) {
    if (eLow <= 0.0 || eHigh <= 0.0) return 0.0;
    return std::sqrt(eLow * eHigh);
}

double PostProcessing::EffAreaErrFromCounts(double n0, double n, double effArea) {
    if (n0 <= 0.0 || n <= 0.0) return 0.0;
    double val = (n0 - n) / (n * n0);
    if (val <= 0.0) return 0.0;
    return effArea * std::sqrt(val);
}

static short GetColorForParticle(const std::string& particleName) {
    if (particleName == "gamma") return kGreen + 2;
    if (particleName == "e-") return kOrange + 7;
    if (particleName == "proton") return kBlue + 1;
    if (particleName == "alpha") return kRed + 1;
    return kBlack;
}

void PostProcessing::SaveHistPng(const std::string& histName,
                                 const std::string& outPngPath,
                                 const std::string& plotTitle,
                                 const std::string& yTitle,
                                 const bool filled,
                                 const bool useTrig) {
    TH1* h = GetHistOrThrow(histName);

    TCanvas canvas("canvas", "canvas", 1920, 1080);
    canvas.SetLogx(true);

    const short color = GetColorForParticle(particleName);

    h->SetTitle(plotTitle.c_str());
    h->GetXaxis()->SetTitle("Energy");
    h->GetYaxis()->SetTitle(yTitle.c_str());

    // h->GetXaxis()->SetMoreLogLabels(true);
    // h->GetXaxis()->SetNoExponent(true);

    if (eMinMeV > 0.0 && eMaxMeV > eMinMeV) {
        h->GetXaxis()->SetRangeUser(eMinMeV, eMaxMeV);
        if (useTrig && eTrigMeV > 0.0)
            h->GetXaxis()->SetRangeUser(eTrigMeV, eMaxMeV);
    }

    h->SetLineColor(color);
    h->SetLineWidth(filled ? 1 : 2);

    if (filled) {
        h->SetFillColorAlpha(color, 0.35);
        h->SetFillStyle(1001);

        h->SetMarkerColor(color);
    } else {
        h->SetFillStyle(0);
    }

    h->Draw("HIST");

    canvas.SaveAs(outPngPath.c_str());
}


void PostProcessing::ExportTreeToCsv(const std::string& treeName,
                                     const std::string& csvPath) {
    TTree* tree = nullptr;
    rootFile->GetObject(treeName.c_str(), tree);
    if (!tree) {
        throw std::runtime_error("TTree/NTuple not found: " + treeName);
    }

    auto* leaves = tree->GetListOfLeaves();
    if (!leaves || leaves->GetEntries() == 0) {
        throw std::runtime_error("No leaves found in tree: " + treeName);
    }

    std::ofstream out(csvPath);
    if (!out.is_open()) {
        throw std::runtime_error("Cannot open output CSV: " + csvPath);
    }


    const int nLeaves = leaves->GetEntries();
    for (int i = 0; i < nLeaves; ++i) {
        auto* leaf = dynamic_cast<TLeaf*>(leaves->At(i));
        out << (leaf ? leaf->GetName() : "unknown");
        if (i + 1 != nLeaves) out << ",";
    }
    out << "\n";


    const Long64_t nEntries = tree->GetEntries();
    for (Long64_t entry = 0; entry < nEntries; ++entry) {
        tree->GetEntry(entry);

        for (int i = 0; i < nLeaves; ++i) {
            auto* leaf = dynamic_cast<TLeaf*>(leaves->At(i));
            if (!leaf) {
                out << "";
                if (i + 1 != nLeaves) out << ",";
                continue;
            }

            if (leaf->InheritsFrom(TLeafC::Class())) {
                auto* lc = dynamic_cast<TLeafC*>(leaf);
                std::string s = lc->GetValueString();

                bool needQuotes = s.find(',') != std::string::npos ||
                    s.find('"') != std::string::npos ||
                    s.find('\n') != std::string::npos;
                if (needQuotes) {
                    std::string esc;
                    esc.reserve(s.size() + 8);
                    for (char c : s) {
                        if (c == '"') esc += "\"\"";
                        else esc += c;
                    }
                    out << "\"" << esc << "\"";
                } else {
                    out << s;
                }
            } else {
                int nData = leaf->GetNdata();
                if (nData <= 1) {
                    double v = leaf->GetValue(0);
                    out << std::setprecision(17) << v;
                } else {
                    std::ostringstream cell;
                    cell << std::setprecision(17);
                    for (int k = 0; k < nData; ++k) {
                        if (k) cell << ";";
                        cell << leaf->GetValue(k);
                    }
                    out << "\"" << cell.str() << "\"";
                }
            }

            if (i + 1 != nLeaves) out << ",";
        }
        out << "\n";
    }

    out.close();
}

void PostProcessing::ExtractNtData() {
    ExportTreeToCsv("edep", (fs::path(csvDir) / "edep.csv").string());
    ExportTreeToCsv("primary", (fs::path(csvDir) / "primary.csv").string());

    if (processSecondaries) {
        ExportTreeToCsv("event", (fs::path(csvDir) / "event.csv").string());
        ExportTreeToCsv("interactions", (fs::path(csvDir) / "interactions.csv").string());
    }

    if (processOptics) {
        ExportTreeToCsv("sipm_event", (fs::path(csvDir) / "sipm_event.csv").string());
        ExportTreeToCsv("sipm_ch", (fs::path(csvDir) / "sipm_ch.csv").string());
    }
}

void PostProcessing::SaveEffArea() {
    TH1* gen = GetHistOrThrow("genEnergyHist");
    TH1* trig = GetHistOrThrow("trigEnergyHist");
    TH1* effArea = GetHistOrThrow("effAreaHist");

    int nBins = gen->GetNbinsX();
    if (trig->GetNbinsX() != nBins || effArea->GetNbinsX() != nBins) {
        throw std::runtime_error("Histogram binning mismatch among genEnergyHist/trigEnergyHist/effAreaHist");
    }


    SaveHistPng("effAreaHist", (fs::path(effectiveAreaDir) / "effective_area.png").string(),
                "Effective Area vs Energy", "Effective Area [cm^{2}]", false, true);

    SaveHistPng("genEnergyHist", (fs::path(histogramsDir) / "genEnergyHist.png").string(),
                "Initial Energy", "Counts", true, false);

    SaveHistPng("trigEnergyHist", (fs::path(histogramsDir) / "trigEnergyHist.png").string(),
                "N_{trig} vs Energy", "Counts", true, false);


    std::ofstream out((fs::path(effectiveAreaDir) / "effective_area_by_energy.csv").string());
    if (!out.is_open()) {
        throw std::runtime_error("Cannot open output CSV: effective_area_by_energy.csv");
    }

    out << "E_low,E_high,E_width,E_center_geom,N0_i,N_i,effective_area,effective_area_err\n";
    out << std::setprecision(17);

    auto* ax = gen->GetXaxis();

    for (int i = 1; i <= nBins; ++i) {
        double eLow = ax->GetBinLowEdge(i);
        double eHigh = ax->GetBinUpEdge(i);
        double eWidth = eHigh - eLow;
        double eCenterGeom = GeomCenter(eLow, eHigh);

        double n0 = gen->GetBinContent(i);
        double n = trig->GetBinContent(i);

        double aeff = effArea->GetBinContent(i);
        double aeffErr = EffAreaErrFromCounts(n0, n, aeff);

        out << eLow << ","
            << eHigh << ","
            << eWidth << ","
            << eCenterGeom << ","
            << n0 << ","
            << n << ","
            << aeff << ","
            << aeffErr << "\n";
    }

    out.close();
}

void PostProcessing::SaveSensitivity() {
    TH1* gen = GetHistOrThrow("genEnergyHist");
    TH1* trig = GetHistOrThrow("trigEnergyHist");
    TH1* sens = GetHistOrThrow("sensitivityHist");

    int nBins = gen->GetNbinsX();
    if (trig->GetNbinsX() != nBins || sens->GetNbinsX() != nBins) {
        throw std::runtime_error("Histogram binning mismatch among genEnergyHist/trigEnergyHist/sensitivityHist");
    }

    SaveHistPng("sensitivityHist", (fs::path(sensitivityDir) / "sensitivity.png").string(),
                "Sensitivity vs Energy", "Sensitivity [cm^{2} #cdot sr]", false, true);

    std::ofstream out((fs::path(sensitivityDir) / "sensitivity_by_energy.csv").string());
    if (!out.is_open()) {
        throw std::runtime_error("Cannot open output CSV: sensitivity_by_energy.csv");
    }

    out << "E_low,E_high,E_width,E_center_geom,N0_i,N_i,sensitivity\n";
    out << std::setprecision(17);

    auto* ax = gen->GetXaxis();

    for (int i = 1; i <= nBins; ++i) {
        double eLow = ax->GetBinLowEdge(i);
        double eHigh = ax->GetBinUpEdge(i);
        double eWidth = eHigh - eLow;
        double eCenterGeom = GeomCenter(eLow, eHigh);

        double n0 = gen->GetBinContent(i);
        double n = trig->GetBinContent(i);

        double s = sens->GetBinContent(i);

        out << eLow << ","
            << eHigh << ","
            << eWidth << ","
            << eCenterGeom << ","
            << n0 << ","
            << n << ","
            << s << "\n";
    }

    out.close();
}

void PostProcessing::SaveTrigEdepCsv() {
    TTree* primary = nullptr;
    rootFile->GetObject("primary", primary);
    if (!primary) {
        throw std::runtime_error("TTree not found: primary");
    }

    Int_t eventID_p = 0;
    double E0 = 0.0;

    primary->SetBranchStatus("*", false);
    primary->SetBranchStatus("eventID", true);
    primary->SetBranchStatus("E_MeV", true);

    primary->SetBranchAddress("eventID", &eventID_p);
    primary->SetBranchAddress("E_MeV", &E0);

    std::unordered_map<int, double> e0ByEvent;
    e0ByEvent.reserve(std::max<Long64_t>(1, primary->GetEntries()));

    const Long64_t nP = primary->GetEntries();
    for (Long64_t i = 0; i < nP; ++i) {
        primary->GetEntry(i);
        e0ByEvent[eventID_p] = E0;
    }

    TTree* edep = nullptr;
    rootFile->GetObject("edep", edep);
    if (!edep) {
        throw std::runtime_error("TTree not found: edep");
    }

    Int_t eventID_e = 0;

    Char_t det_name[64] = {};

    double edep_MeV = 0.0;

    edep->SetBranchStatus("*", false);
    edep->SetBranchStatus("eventID", true);
    edep->SetBranchStatus("det_name", true);
    edep->SetBranchStatus("edep_MeV", true);

    edep->SetBranchAddress("eventID", &eventID_e);
    edep->SetBranchAddress("det_name", det_name);
    edep->SetBranchAddress("edep_MeV", &edep_MeV);

    struct Agg {
        double crystal = 0.0;
        double veto = 0.0;
        double bottomVeto = 0.0;
    };

    std::unordered_map<int, Agg> agg;
    agg.reserve(std::max<Long64_t>(1, edep->GetEntries()));

    const Long64_t nE = edep->GetEntries();
    for (Long64_t i = 0; i < nE; ++i) {
        edep->GetEntry(i);

        auto& a = agg[eventID_e];

        if (std::strcmp(det_name, "Crystal") == 0) {
            a.crystal += edep_MeV;
        } else if (std::strcmp(det_name, "Veto") == 0) {
            a.veto += edep_MeV;
        } else if (std::strcmp(det_name, "BottomVeto") == 0) {
            a.bottomVeto += edep_MeV;
        }
    }

    const std::string outPath = (fs::path(histogramsDir) / "trig_edep.csv").string();
    std::ofstream out(outPath);
    if (!out.is_open()) {
        throw std::runtime_error("Cannot open output CSV: " + outPath);
    }

    out << "eventID,E0,Crystal_only_edep\n";
    out << std::setprecision(17);

    std::vector<int> events;
    events.reserve(e0ByEvent.size());
    for (const auto& kv : e0ByEvent) events.push_back(kv.first);
    std::sort(events.begin(), events.end());

    for (int evt : events) {
        const double e0 = e0ByEvent.at(evt);

        double crystal_only = 0.0;
        auto it = agg.find(evt);
        if (it != agg.end()) {
            const double c = it->second.crystal;
            const double v = it->second.veto;
            const double b = it->second.bottomVeto;

            if (c > 0.0 && v + b == 0.0) {
                crystal_only = c;
            }
        }

        out << evt << "," << e0 << "," << crystal_only << "\n";
    }

    out.close();
}
