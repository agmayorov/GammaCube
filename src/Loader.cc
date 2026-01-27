#include "Loader.hh"

Loader::Loader(int argc, char** argv) {
    numThreads = G4Threading::G4GetNumberOfCores();
    useUI = true;
    macroFile = "../run.mac";
    detectorType = "CsI";
    fluxType = "Uniform";
    geomConfigPath = "../geometry_config.txt";
    fluxDirection = "isotropic";
    useOptics = false;
    viewDeg = 360 * deg;
    yieldScale = 1;
    nBins = 1000;
    outputFile = "GammaCube.root";
    eCrystalThreshold = 0 * MeV;
    eVetoThreshold = 0 * MeV;
    saveSecondaries = false;

    for (int i = 0; i < argc; i++) {
        if (std::string input(argv[i]); input == "-i" || input == "-input") {
            macroFile = argv[i + 1];
            useUI = false;
            viewDeg = 360 * deg;
        } else if (input == "-t" || input == "-threads") {
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
        } else if (input == "-noUI") {
            useUI = false;
        } else if (input == "-d" || input == "--detector") {
            detectorType = argv[i + 1];
        } else if (input == "-f" || input == "--flux-type") {
            fluxType = argv[i + 1];
        } else if (input == "--flux-dir" || input == "--f-dir" || input == "-fd") {
            fluxDirection = argv[i + 1];
        } else if (input == "--use-optics") {
            useOptics = true;
        } else if (input == "--save-secondaries") {
            saveSecondaries = true;
        } else if (input == "-g" || input == "--geom-config") {
            geomConfigPath = argv[i + 1];
        } else if (input == "-o" || input == "--output-file") {
            outputFile = argv[i + 1];
            outputFile += ".root";
        }
    }

    configPath = "../Flux_config/" + fluxType + "_params.txt";

    CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);
    CLHEP::HepRandom::setTheSeed(time(nullptr));

#ifdef G4MULTITHREADED
    runManager = new G4MTRunManager;
    runManager->SetNumberOfThreads(numThreads);
#else
    runManager = new G4RunManager;
#endif

    auto* realWorld = new Geometry(detectorType, useOptics, viewDeg, yieldScale);
    runManager->SetUserInitialization(realWorld);
    auto* physicsList = new FTFP_BERT;
    physicsList->ReplacePhysics(new G4EmStandardPhysics_option4());
    physicsList->ReplacePhysics(new G4RadioactiveDecayPhysics());

    if (useOptics) {
        auto* opticalPhysics = new G4OpticalPhysics();

        auto* op = G4OpticalParameters::Instance();
        op->SetProcessActivation("Cerenkov", false);
        op->SetProcessActivation("Scintillation", true);
        op->SetProcessActivation("OpAbsorption", true);
        op->SetProcessActivation("OpRayleigh", true);
        op->SetProcessActivation("OpMieHG", false);
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

    G4double EminMeV = std::max({std::stod(ReadValue("E_min:", "")) * MeV, eCrystalThreshold});
    G4double EmaxMeV = std::stod(ReadValue("E_max:", "")) * MeV;

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
    area = Area_cm2(Sizes::modelRadius, Sizes::modelHeight, dir);
    runManager->SetUserInitialization(new ActionInitialization(fluxDirection, fluxType, useOptics, saveSecondaries,
                                                               eCrystalThreshold, eVetoThreshold, area, nBins, EminMeV,
                                                               EmaxMeV, outputFile));
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

    if (fluxType == "PLAW") {
        fType = FluxType::PLAW;
        fp.A = std::stod(ReadValue("A:"));
        fp.alpha = std::stod(ReadValue("alpha:"));
        fp.E_piv = std::stod(ReadValue("E_Piv:"));
        er.Emin = std::stod(ReadValue("E_min:"));
        er.Emax = std::stod(ReadValue("E_max:"));
    } else if (fluxType == "COMP") {
        fType = FluxType::COMP;
        fp.A = std::stod(ReadValue("A:"));
        fp.alpha = std::stod(ReadValue("alpha:"));
        fp.E_piv = std::stod(ReadValue("E_Piv:"));
        fp.E_peak = std::stod(ReadValue("E_Peak:"));
        er.Emin = std::stod(ReadValue("E_min:"));
        er.Emax = std::stod(ReadValue("E_max:"));
    } else if (fluxType == "SEP") {
        fType = FluxType::SEP;
        fp.sep_year = std::stoi(ReadValue("year:"));
        fp.sep_order = std::stoi(ReadValue("order:"));
        fp.sep_csv_path = "../SEP_coefficients.CSV";
        er.Emin = std::stod(ReadValue("E_min:"));
        er.Emax = std::stod(ReadValue("E_max:"));
    } else if (fluxType == "Galactic") {
        fType = FluxType::GALACTIC;
        fp.phiMV = std::stod(ReadValue("phiMV:"));
        fp.particle = ReadValue("particle:");
        er.Emin = std::stod(ReadValue("E_min:"));
        er.Emax = std::stod(ReadValue("E_max:"));
    } else if (fluxType == "Table") {
        fType = FluxType::TABLE;
        fp.particle = ReadValue("particle:");
        fp.table_path = ReadValue("table_path:");
        er.Emin = std::stod(ReadValue("E_min:"));
        er.Emax = std::stod(ReadValue("E_max:"));
    } else {
        fType = FluxType::UNIFORM;
        er.Emin = std::stod(ReadValue("E_min:"));
        er.Emax = std::stod(ReadValue("E_max:"));
    }

    RateCounts counts{crystalOnly, crystalAndVeto};

    RateResult rr{};
    bool rate_ok = true;
    try {
        rr = computeRate(fType, fp, er, area, N, counts);
    }
    catch (const std::exception& ex) {
        rate_ok = false;
        G4cerr << "[SaveConfig] WARNING: rate computation failed: " << ex.what() << G4endl;
    }

    RateResult rrReal{};
    bool rate_real_ok = true;
    try {
        rrReal = computeRateReal(fType, fp, er, effArea, nBins);
    }
    catch (const std::exception& ex) {
        rate_real_ok = false;
        G4cerr << "[SaveConfig] WARNING: Real rate computation failed: " << ex.what() << G4endl;
    }

    std::ostringstream buf;

    buf << "N: " << N << "\n\n";
    buf << "Detector_type: " << detectorType << "\n\n";
    buf << "Flux_type: " << fluxType << "\n";
    buf << "Flux_dir: " << fluxDirection << "\n";

    buf << "Flux_params:\n{\n\t";
    if (fluxType == "PLAW") {
        buf << "A: " << std::stod(ReadValue("A:")) << ",\n\t";
        buf << "alpha: " << std::stod(ReadValue("alpha:")) << ",\n\t";
        buf << "E_Piv: " << std::stod(ReadValue("E_Piv:")) << "\n";
    } else if (fluxType == "COMP") {
        buf << "A: " << std::stod(ReadValue("A:")) << ",\n\t";
        buf << "alpha: " << std::stod(ReadValue("alpha:")) << ",\n\t";
        buf << "E_Piv: " << std::stod(ReadValue("E_Piv:")) << ",\n\t";
        buf << "E_Peak: " << std::stod(ReadValue("E_Peak:")) << "\n";
    } else if (fluxType == "SEP") {
        buf << "year: " << std::stoi(ReadValue("year:")) << ",\n\t";
        buf << "order: " << std::stoi(ReadValue("order:")) << "\n";
    } else if (fluxType == "Galactic") {
        buf << "phiMV: " << std::stod(ReadValue("phiMV:")) << ",\n\t";
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
            buf << (i == 0 ? "" : "\t") << particles[i] << ": (" << EminVec[i] << ", " << EmaxVec[i] << "),\n";
        }
    }
    if (fluxType != "Uniform") {
        buf << "(" << ReadValue("E_min:") << ", " << ReadValue("E_max:") << ")\n";
    }
    buf << "}\n\n";

    buf << "Counts:\n{\n\t";
    buf << "Crystal_only: " << crystalOnly << "\n\t";
    buf << "Veto_then_Crystal: " << crystalAndVeto << "\n}\n\n";

    buf << "Thresholds:\n{\n\t";
    buf << std::fixed << std::setprecision(6);
    if (rate_ok) {
        buf << "Crystal: " << eCrystalThreshold << "\n\t";
        buf << "Veto: " << eVetoThreshold << "\n";
    }
    buf << "}\n\n";

    buf << "Rates:\n{\n\t";
    buf << std::fixed << std::setprecision(6);
    if (rate_ok) {
        buf << "Area: " << area << "\n\t";
        buf << "Integral: " << rr.integral << "\n\t";
        buf << "Ndot: " << rr.Ndot << "\n\t";
        buf << "Rate_Crystal_only: " << rr.rateCrystal << "\n\t";
        buf << "Rate_Both: " << rr.rateBoth << "\n\t";
    } else {
        buf << "Area: NaN\n\t";
        buf << "Integral: NaN\n\t";
        buf << "Ndot: NaN\n\t";
        buf << "Rate_Crystal_only: NaN\n\t";
        buf << "Rate_Both: NaN\n";
    }
    if (rate_real_ok) {
        buf << "Rate_Real: " << rrReal.rateRealCrystal << "\n";
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

    G4double Emin = std::max({std::stod(ReadValue("E_min:")), eCrystalThreshold});
    G4double Emax = std::stod(ReadValue("E_max:"));
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
        PostProcessing postProcessing(outputFile, outDir, saveSecondaries, useOptics, Emin, Emax, eCrystalThreshold,
                                      part);

        postProcessing.ExtractNtData();
        if (Emin < Emax) {
            if (fluxDirection == "isotropic" || fluxDirection == "isotropic_down" || fluxDirection == "isotropic_up")
                postProcessing.SaveSensitivity();
            else
                postProcessing.SaveEffArea();
        }
        postProcessing.SaveTrigEdepCsv();
        postProcessing.SaveEdepCsv();

        std::cout << "Done!\n";
    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
    }
}
