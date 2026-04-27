#include "Loader.hh"

using namespace Configuration;

Loader::Loader(int argc, char** argv) {
    numThreads = G4Threading::G4GetNumberOfCores();
    useUI = true;
    macroFile = "../run.mac";
    geomConfigPath = "../geometry_txt";
    detectorType = "CsI";
    fluxType = "Uniform";
    fluxDirection = "isotropic";
    eCrystalThreshold = 0 * MeV;
    eVetoThreshold = 0 * MeV;
    useOptics = false;
    yieldScale = 1;
    oCrystalThreshold = 0 * MeV;
    oVetoThreshold = 0 * MeV;
    outputFile = "GammaCube.root";
    nBins = 1000;
    saveSecondaries = false;
    savePhotons = false;

    for (int i = 0; i < argc; i++) {
        if (std::string input(argv[i]); input == "-i" || input == "--input") {
            macroFile = argv[i + 1];
            useUI = false;
            viewDeg = 360 * deg;
        } else if (input == "-t" || input == "--threads") {
            numThreads = std::stoi(argv[i + 1]);
        } else if (input == "-ys" || input == "--yield-scale") {
            yieldScale = std::stoi(argv[i + 1]);
        } else if (input == "--bins") {
            nBins = std::stoi(argv[i + 1]);
        } else if ((input == "-vd" || input == "--view-deg") and useUI) {
            viewDeg = std::stod(argv[i + 1]) * deg;
        } else if (input == "-ct" || input == "--crystal-threshold") {
            eCrystalThreshold = std::stod(argv[i + 1]) * MeV;
        } else if (input == "-vt" || input == "--veto-threshold") {
            eVetoThreshold = std::stod(argv[i + 1]) * MeV;
        } else if (input == "-oct" || input == "--crystal-optic-threshold") {
            oCrystalThreshold = std::stoi(argv[i + 1]);
        } else if (input == "-ovt" || input == "--veto-optic-threshold") {
            oVetoThreshold = std::stoi(argv[i + 1]);
        } else if (input == "-obvt" || input == "--bottom-veto-optic-threshold") {
            oBottomVetoThreshold = std::stoi(argv[i + 1]);
        } else if ((input == "-vch" || input == "--veto-chamfer-height") and Sizes::vetoTopRoundedRadius == 0) {
            Sizes::vetoChamferHeight = std::stod(argv[i + 1]) * mm;
        } else if ((input == "-vtr" || input == "--veto-top-rounded") and Sizes::vetoTopRoundedRadius == 0) {
            Sizes::vetoTopRoundedRadius = std::stod(argv[i + 1]) * mm;
            if (Sizes::vetoTopRoundedRadius > 0)
                Sizes::vetoChamferHeight = 0 * mm;
        } else if (input == "-noUI") {
            useUI = false;
        } else if (input == "--polished") {
            polishedTyvek = true;
        } else if (input == "-d" || input == "--detector") {
            detectorType = argv[i + 1];
        } else if (input == "-csc" || input == "--crystal-sipm-config" || input == "-sipm") {
            crystalSiPMConfig = argv[i + 1];
        } else if (input == "-f" || input == "--flux-type") {
            fluxType = argv[i + 1];
        } else if (input == "--flux-dir" || input == "--f-dir" || input == "-fd") {
            fluxDirection = argv[i + 1];
        } else if (input == "--use-optics") {
            useOptics = true;
        } else if (input == "--save-secondaries") {
            saveSecondaries = true;
        } else if (input == "--save-photons") {
            savePhotons = true;
        } else if (input == "-g" || input == "--geom-config") {
            geomConfigPath = argv[i + 1];
        } else if (input == "-w4" || input == "--weight4") {
            weight4 = std::stod(argv[i + 1]);
        } else if (input == "-w8" || input == "--weight8") {
            weight8 = std::stod(argv[i + 1]);
        } else if (input == "-o" || input == "--output-file") {
            outputFile = argv[i + 1];
            outputFile += ".root";
        }
    }

    savePhotons = savePhotons and useOptics;

    configPath = "../Flux_config/" + fluxType + "_params.txt";

    CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);
    CLHEP::HepRandom::setTheSeed(time(nullptr));

#ifdef G4MULTITHREADED
    runManager = new G4MTRunManager;
    runManager->SetNumberOfThreads(numThreads);
#else
    runManager = new G4RunManager;
#endif

    auto* realWorld = new Geometry();
    runManager->SetUserInitialization(realWorld);
    auto* physicsList = new FTFP_BERT;
    physicsList->ReplacePhysics(new G4EmStandardPhysics_option4());
    physicsList->ReplacePhysics(new G4RadioactiveDecayPhysics());

    if (useOptics) {
        auto* opticalPhysics = new G4OpticalPhysics();

        auto* op = G4OpticalParameters::Instance();
        op->SetProcessActivation("Cerenkov", true);
        op->SetProcessActivation("Scintillation", true);
        op->SetProcessActivation("OpAbsorption", true);
        op->SetProcessActivation("OpRayleigh", true);
        op->SetProcessActivation("OpMieHG", true);
        op->SetProcessActivation("OpBoundary", true);

        op->SetScintTrackSecondariesFirst(true);
        op->SetCerenkovTrackSecondariesFirst(false);
        // op->SetScintByParticleType(true);
        // op->SetCerenkovMaxPhotonsPerStep(200);
        // op->SetCerenkovMaxBetaChange(10.0);

        physicsList->RegisterPhysics(opticalPhysics);
    }

    physicsList->RegisterPhysics(new G4StepLimiterPhysics());
    runManager->SetUserInitialization(physicsList);

    Emin = std::max({std::stod(ReadValue("E_min:", "")) * MeV, eCrystalThreshold});
    Emax = std::stod(ReadValue("E_max:", "")) * MeV;
    if (fluxType == "Table") {
        std::string path = ReadValue("table_path:", "");
        auto energyTable = Utils::ReadCSV(path, 1., false, MeV);

        Emin = std::max({energyTable.GetMinE(), Emin});
        Emax = std::min({energyTable.GetMaxE(), Emax});
    }

    if (fluxDirection == "isotropic") {
        dir = FluxDir::Isotropic;
    } else if (fluxDirection == "isotropic_up") {
        dir = FluxDir::Isotropic_up;
    } else if (fluxDirection == "isotropic_down") {
        dir = FluxDir::Isotropic_down;
    } else if (fluxDirection == "vertical_up") {
        dir = FluxDir::Vertical_up;
    } else if (fluxDirection == "vertical_down") {
        dir = FluxDir::Vertical_down;
    } else if (fluxDirection == "horizontal") {
        dir = FluxDir::Horizontal;
    }
    if (fluxType == "Uniform") {
        std::string isLogStr = ReadValue("is_log:", configPath);
        isLogBin = isLogStr == "1" || isLogStr == "true";
    }
    std::cout << isLogBin << std::endl;
    area = Area_cm2(Sizes::modelRadius, Sizes::modelHeight, dir);
    runManager->SetUserInitialization(new ActionInitialization(area));
    runManager->Initialize();

    visManager = new G4VisExecutive;
    visManager->Initialize();
    G4UImanager* UImanager = G4UImanager::GetUIpointer();

    if (!useUI) {
        const G4String command = "/control/execute ";
        UImanager->ApplyCommand(command + macroFile);
    } else {
        auto* ui = new G4UIExecutive(argc, argv, "qt");
        UImanager->ApplyCommand("/control/execute ../vis.mac");
        ui->SessionStart();
        delete ui;
    }

    const auto* runAction = dynamic_cast<const RunAction*>(runManager->GetUserRunAction());
    if (runAction) {
        const auto& [cOnly, cAndV] = runAction->GetCounts();
        crystalOnly = cOnly;
        crystalAndVeto = cAndV;
        effArea = runAction->GetEffArea();
        const auto& [cOnlyOpt, cAndVOpt] = runAction->GetOptCounts();
        crystalOnlyOpt = cOnlyOpt;
        crystalAndVetoOpt = cAndVOpt;
        effAreaOpt = runAction->GetEffAreaOpt();
    }
    SaveConfig();
    RunPostProcessing();
}

Loader::~Loader() {
    delete runManager;
    delete visManager;
}


std::string Loader::ReadValue(const std::string& key, const std::string& filepath = "") const {
    std::ifstream file(filepath.empty() ? configPath : filepath);
    if (!file.is_open()) {
        G4Exception("Loader::ReadValue", "FILE_OPEN_FAIL",
                    FatalException, ("Cannot open " + (filepath.empty() ? configPath : filepath)).c_str());
    }

    std::string line;
    while (std::getline(file, line)) {
        if (line.find(key) != std::string::npos) {
            return line.substr(key.length() + 1);
        }
    }

    return "";
}


inline std::string Trim(std::string st) {
    auto notSpace = [](const unsigned char c) {
        return !std::isspace(c);
    };
    st.erase(st.begin(), std::find_if(st.begin(), st.end(), notSpace));
    st.erase(std::find_if(st.rbegin(), st.rend(), notSpace).base(), st.end());
    return st;
}


std::vector<G4String> Split(const G4String& line) {
    std::vector<G4String> result;
    std::stringstream ss(line);
    G4String token;
    while (std::getline(ss, token, ',')) {
        token = Trim(token);
        if (!token.empty())
            result.push_back(token);
    }
    return result;
}


void Loader::SaveConfig() const {
    const int N = std::stoi(ReadValue("/run/beamOn", "../run.mac"));

    EnergyRange er{};
    FluxType fType{};
    FluxParams fp{};

    er.Emin = Emin;
    er.Emax = Emax;

    if (fluxType == "PLAW") {
        fType = FluxType::PLAW;
        fp.A = std::stod(ReadValue("A:"));
        fp.alpha = std::stod(ReadValue("alpha:"));
        fp.E_piv = std::stod(ReadValue("E_Piv:"));
    } else if (fluxType == "COMP") {
        fType = FluxType::COMP;
        fp.A = std::stod(ReadValue("A:"));
        fp.alpha = std::stod(ReadValue("alpha:"));
        fp.E_piv = std::stod(ReadValue("E_Piv:"));
        fp.E_peak = std::stod(ReadValue("E_Peak:"));
    } else if (fluxType == "SEP") {
        fType = FluxType::SEP;
        fp.sep_year = std::stoi(ReadValue("year:"));
        fp.sep_order = std::stoi(ReadValue("order:"));
        fp.sep_csv_path = "../SEP_coefficients.CSV";
    } else if (fluxType == "Galactic") {
        fType = FluxType::GALACTIC;
        fp.phiMV = std::stod(ReadValue("phiMV:"));
        fp.particle = ReadValue("particle:");
    } else if (fluxType == "Table") {
        fType = FluxType::TABLE;
        fp.particle = ReadValue("particle:");
        fp.table_path = ReadValue("table_path:");
    } else {
        fType = FluxType::UNIFORM;
    }

    RateCounts counts{crystalOnly, crystalAndVeto};
    RateCounts countsOpt{crystalOnlyOpt, crystalAndVetoOpt};

    RateResult rr{};
    bool rate_ok = true;
    try {
        rr = computeRate(fType, fp, er, area, N, counts);
    }
    catch (const std::exception& ex) {
        rate_ok = false;
    }

    RateResult rrReal{};
    bool rate_real_ok = true;
    try {
        rrReal = computeRateReal(fType, fp, er, effArea, nBins);
    }
    catch (const std::exception& ex) {
        rate_real_ok = false;
    }

    RateResult rr_opt{};
    bool rate_opt_ok = useOptics;
    try {
        rr_opt = computeRate(fType, fp, er, area, N, countsOpt);
    }
    catch (const std::exception& ex) {
        rate_opt_ok = false;
        // G4cerr << "[SaveConfig] WARNING: rate computation failed: " << ex.what() << G4endl;
    }

    RateResult rrReal_opt{};
    bool rate_real_opt_ok = useOptics;
    try {
        rrReal_opt = computeRateReal(fType, fp, er, effAreaOpt, nBins);
    }
    catch (const std::exception& ex) {
        rate_real_opt_ok = false;
        // G4cerr << "[SaveConfig] WARNING: Real rate computation failed: " << ex.what() << G4endl;
    }

    std::ostringstream buf;

    buf << "N: " << N << "\n\n";
    buf << "Detector_type: " << detectorType << "\n";
    buf << "Crystal_SiPM_configuration: " << crystalSiPMConfig << "\n";
    buf << "Tyvek_surface: " << (polishedTyvek ? "polished" : "diffuse") << "\n\n";
    buf << "Use_optics: " << useOptics << "\n\n";
    buf << "Flux_type: " << fluxType << "\n";
    buf << "Flux_dir: " << fluxDirection << "\n";

    buf << "Flux_params:\n{\n\t";
    if (fluxType == "PLAW") {
        buf << "A: " << std::stod(ReadValue("A:")) << ",\n\t";
        buf << "alpha: " << std::stod(ReadValue("alpha:")) << ",\n\t";
        buf << "E_Piv: " << std::stod(ReadValue("E_Piv:")) << " MeV\n";
    } else if (fluxType == "COMP") {
        buf << "A: " << std::stod(ReadValue("A:")) << ",\n\t";
        buf << "alpha: " << std::stod(ReadValue("alpha:")) << ",\n\t";
        buf << "E_Piv: " << std::stod(ReadValue("E_Piv:")) << " MeV,\n\t";
        buf << "E_Peak: " << std::stod(ReadValue("E_Peak:")) << " MeV\n";
    } else if (fluxType == "SEP") {
        buf << "year: " << std::stoi(ReadValue("year:")) << ",\n\t";
        buf << "order: " << std::stoi(ReadValue("order:")) << "\n";
    } else if (fluxType == "Galactic") {
        buf << "phiMV: " << std::stod(ReadValue("phiMV:")) << " MV,\n\t";
        buf << "particle: " << ReadValue("particle:") << "\n";
    } else if (fluxType == "Table") {
        buf << "table_path: " << ReadValue("table_path:") << ",\n\t";
        buf << "particle: " << ReadValue("particle:") << "\n";
    } else if (fluxType == "Uniform") {
        buf << "fractions: " << ReadValue("fractions:") << "\n";
    }
    buf << "}\n\n";

    buf << "Particles: [";
    if (fluxType == "PLAW") {
        buf << "gamma";
    } else if (fluxType == "SEP") {
        buf << "proton";
    } else if (fluxType == "Galactic" or fluxType == "Table") {
        buf << ReadValue("particle:");
    } else if (fluxType == "Uniform") {
        buf << ReadValue("particles:");
    }
    buf << "]\n";

    buf << "Energies:\n{\n\t";
    if (fluxType == "PLAW" || fluxType == "COMP") {
        buf << "gamma: ";
    } else if (fluxType == "SEP") {
        buf << "proton: ";
    } else if (fluxType == "Galactic" or fluxType == "Table") {
        buf << ReadValue("particle:") << ": ";
    } else if (fluxType == "Uniform") {
        std::vector<G4String> particles = Split(ReadValue("particles:"));
        std::vector<G4String> EminVec = Split(ReadValue("E_min:"));
        std::vector<G4String> EmaxVec = Split(ReadValue("E_max:"));
        for (size_t i = 0; i < particles.size(); i++) {
            buf << (i == 0 ? "" : "\t") << particles[i] << ": (" << EminVec[i] << " MeV, " << EmaxVec[i] << " MeV),\n";
        }
    }
    if (fluxType != "Uniform") {
        buf << "(" << Emin << " MeV, " << Emax << " MeV)\n";
    }
    buf << "}\n\n";

    buf << "Counts:\n{\n\t";
    buf << "Crystal_only: " << crystalOnly << "\n\t";
    buf << "Veto_then_Crystal: " << crystalAndVeto << "\n}\n\n";

    buf << "Optical_Counts:\n{\n\t";
    buf << "Crystal_only: " << crystalOnlyOpt << "\n\t";
    buf << "Veto_then_Crystal: " << crystalAndVetoOpt << "\n}\n\n";

    buf << "Thresholds:\n{\n\t";
    if (rate_ok) {
        buf << "Crystal: " << eCrystalThreshold << " MeV\n\t";
        buf << "Veto: " << eVetoThreshold << " MeV\n";
    }
    buf << "}\n\n";

    buf << "Optical_thresholds:\n{\n\t";
    if (rate_ok) {
        buf << "Crystal: " << oCrystalThreshold << "\n\t";
        buf << "Veto: " << oVetoThreshold << "\n\t";
        buf << "BottomVeto: " << oBottomVetoThreshold << "\n";
    }
    buf << "}\n\n";

    if (crystalSiPMConfig == "12-rhombus") {
        buf << "Weight_for_4: " << weight4 << "\n";
        buf << "Weight_for_8: " << weight8 << "\n\n";
    }

    std::string area_dim = fluxDirection.find("isotropic") != std::string::npos ? " sr * cm^2" : " cm^2";
    std::string area_dim_inv = fluxDirection.find("isotropic") != std::string::npos ? " sr^-1 * cm^-2" : " cm^-2";

    buf << "Rates:\n{\n\t";
    buf << std::fixed << std::setprecision(6);
    if (rate_ok) {
        buf << "Area: " << area << area_dim << "\n\t";
        buf << "Integral: " << rr.integral / (fluxType == "Galactic" ? 10000 : 1) << area_dim_inv << " * s^-1\n\t";
        buf << "Ndot: " << rr.Ndot << " s^-1\n\t";
        buf << "Rate_Crystal_only: " << rr.rateCrystal << " s^-1\n\t";
        buf << "Rate_Both: " << rr.rateBoth << " s^-1\n\t";
    } else {
        buf << "Area: NaN\n\t";
        buf << "Integral: NaN\n\t";
        buf << "Ndot: NaN\n\t";
        buf << "Rate_Crystal_only: NaN\n\t";
        buf << "Rate_Both: NaN\n\t";
    }
    if (rate_real_ok) {
        buf << "Rate_Real: " << rrReal.rateRealCrystal << " s^-1\n";
    } else {
        buf << "Rate_Real: NaN\n";
    }
    buf << "}\n\n";

    buf << "Optical_rates:\n{\n\t";
    buf << std::fixed << std::setprecision(6);
    if (rate_opt_ok) {
        buf << "Area: " << area << area_dim << "\n\t";
        buf << "Integral: " << rr_opt.integral / (fluxType == "Galactic" ? 10000 : 1) << area_dim_inv << " * s^-1\n\t";
        buf << "Ndot: " << rr_opt.Ndot << " s^-1\n\t";
        buf << "Rate_Crystal_only: " << rr_opt.rateCrystal << " s^-1\n\t";
        buf << "Rate_Both: " << rr_opt.rateBoth << " s^-1\n\t";
    } else {
        buf << "Area: NaN\n\t";
        buf << "Integral: NaN\n\t";
        buf << "Ndot: NaN\n\t";
        buf << "Rate_Crystal_only: NaN\n\t";
        buf << "Rate_Both: NaN\n\t";
    }
    if (rate_real_opt_ok) {
        buf << "Rate_Real: " << rrReal_opt.rateRealCrystal << " s^-1\n";
    } else {
        buf << "Rate_Real: NaN\n";
    }
    buf << "}\n\n";

    auto sanitize = [](std::string ss) {
        for (char& c : ss) if (c == ' ') c = '_';
        return ss;
    };

    std::string filename = "info_" + detectorType + "_" + fluxType;
    if (fluxType == "Galactic") {
        const std::string part = ReadValue("particle:");
        const std::string phi = ReadValue("phiMV:");
        filename += "_particle:" + part + "_phiMV:" + phi + ".txt";
    } else if (fluxType == "Uniform") {
        const std::string part = ReadValue("particles:");
        filename += "_particle:" + part + ".txt";
    } else {
        filename += ".txt";
    }
    filename = sanitize(filename);

    std::ofstream out(filename);
    if (!out.is_open()) {
        G4cerr << "Ошибка: не удалось открыть файл " << filename << G4endl;
        return;
    }
    out << buf.str();
    out.close();

    std::cout << "Configuration saved in " << filename << std::endl;
}


void Loader::RunPostProcessing() const {
    auto sanitize = [](std::string ss) {
        for (char& c : ss) if (c == ' ') c = '_';
        return ss;
    };
    std::string part;

    try {
        std::cout << "Processing... ";
        std::string outDir = fluxType;
        if (fluxType == "Galactic") {
            const std::string phi = ReadValue("phiMV:");
            part = ReadValue("particle:");
            outDir += "_particle:" + part + "_phiMV:" + phi;
        } else if (fluxType == "Uniform") {
            part = ReadValue("particles:");
            outDir += "_particles:" + part;
        } else if (fluxType == "PLAW" || fluxType == "COMP") {
            part = "gamma";
        } else if (fluxType == "SEP") {
            part = "proton";
        }
        outDir = sanitize(outDir);
        PostProcessing postProcessing(outDir, part);

        postProcessing.ExtractNtData();
        if (Emin < Emax) {
            if (fluxDirection.find("isotropic") != std::string::npos)
                postProcessing.SaveSensitivity();
            else
                postProcessing.SaveEffArea();
        }
        postProcessing.SaveTrigEdepCsv();
        postProcessing.SaveEdepCsv();
        postProcessing.SavePrimaryCsv();
        if (useOptics) {
            postProcessing.SaveOpticsCsv();
        }

        std::cout << "Done!\n";
    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
    }
}
