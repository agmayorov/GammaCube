#include "Loader.hh"

Loader::Loader(int argc, char **argv) {
    numThreads = G4Threading::G4GetNumberOfCores();
    useUI = true;
    macroFile = "../run.mac";
    temperature = 20;
    detectorType = "NaI";
    fluxType = "Uniform";
    sizes.vetoThick = 0 * mm;
    sizes.shellThick = 0 * mm;
    sizes.gapSize = 0 * mm;
    sizes.tyvekThick = 1 * mm;
    sizes.LEDSize = 5 * mm;
    verticalFlux = false;
    doubleLED = false;

    for (int i = 0; i < argc; i++) {
        if (std::string input(argv[i]); input == "-i" || input == "-input") {
            macroFile = argv[i + 1];
            useUI = false;
        } else if (input == "-t" || input == "-threads") {
            numThreads = std::stoi(argv[i + 1]);
        } else if (input == "-noUI") {
            useUI = false;
        } else if (input == "-temp" || input == "--temperature") {
            temperature = std::stod(argv[i + 1]);
        } else if (input == "-s" || input == "--shellThick") {
            sizes.shellThick = std::stod(argv[i + 1]) * mm;
        } else if (input == "-v" || input == "--veto") {
            sizes.vetoThick = std::stod(argv[i + 1]) * mm;
        } else if (input == "-g" || input == "--gap") {
            sizes.gapSize = std::stod(argv[i + 1]) * mm;
        } else if (input == "-tyvek") {
            sizes.tyvekThick = std::stod(argv[i + 1]) * mm;
        } else if (input == "--double-LED") {
            doubleLED = true;
        } else if (input == "-LED") {
            sizes.LEDSize = std::stod(argv[i + 1]);
        } else if (input == "-d" || input == "--detector") {
            detectorType = argv[i + 1];
        } else if (input == "-f" || input == "--flux-type") {
            fluxType = argv[i + 1];
        } else if (input == "--vertical-flux") {
            verticalFlux = true;
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
    Geometry *realWorld = new Geometry(detectorType, sizes, temperature, doubleLED);
    runManager->SetUserInitialization(realWorld);
    G4VModularPhysicsList *physicsList = new FTFP_BERT;
    physicsList->ReplacePhysics(new G4EmStandardPhysics_option4());
    physicsList->ReplacePhysics(new G4RadioactiveDecayPhysics());
    // auto* opticalPhysics = new G4OpticalPhysics();
    // auto* op = G4OpticalParameters::Instance();
    // op->SetProcessActivation("Scintillation", true);
    // op->SetProcessActivation("OpAbsorption", true);
    // op->SetProcessActivation("OpRayleigh", true);
    // op->SetProcessActivation("OpMieHG", true);
    // op->SetProcessActivation("OpBoundary", true);
    // physicsList->RegisterPhysics(opticalPhysics);
    physicsList->RegisterPhysics(new G4StepLimiterPhysics());
    runManager->SetUserInitialization(physicsList);

    runManager->SetUserInitialization(new ActionInitialization(realWorld->sizes, realWorld->modelSize,
                                                               realWorld->detContainerPos, verticalFlux, fluxType));
    runManager->Initialize();

    visManager = new G4VisExecutive;
    visManager->Initialize();
    G4UImanager *UImanager = G4UImanager::GetUIpointer();

    if (!useUI) {
        const G4String command = "/control/execute ";
        UImanager->ApplyCommand(command + macroFile);
    } else {
        G4UIExecutive *ui = new G4UIExecutive(argc, argv, "qt");
        UImanager->ApplyCommand("/control/execute ../vis.mac");
        ui->SessionStart();
        delete ui;
    }

    const auto *ra = dynamic_cast<const RunAction *>(runManager->GetUserRunAction());
    if (ra) {
        const auto &[cOnly, cAndV] = ra->GetCounts();
        crystalOnly = cOnly;
        crystalAndVeto = cAndV;
    }
    SaveConfig(realWorld);
}

Loader::~Loader() {
    delete runManager;
    delete visManager;
}


std::string Loader::ReadValue(const std::string &key, const std::string &filepath = "") const {
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


void Loader::SaveConfig(const Geometry *geometry) const {
    const int N = std::stoi(ReadValue("/run/beamOn", "../run.mac"));

    const double R_mm = geometry->modelSize.y();
    const double H_mm = geometry->modelSize.z() * 2.0;
    const FluxDir dir = verticalFlux ? FluxDir::Vertical : FluxDir::Isotropic;
    double A_eff_cm2 = effectiveArea_cm2(R_mm, H_mm, dir);

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
    } else {
        fType = FluxType::UNIFORM;
        er.Emin = 0.001;
        er.Emax = 500.0;
    }

    RateCounts counts{crystalOnly, crystalAndVeto};

    RateResult rr{};
    bool rate_ok = true;
    try {
        rr = computeRate(fType, fp, er, A_eff_cm2, N, counts);
    } catch (const std::exception &ex) {
        rate_ok = false;
        G4cerr << "[SaveConfig] WARNING: rate computation failed: " << ex.what() << G4endl;
    }

    std::ostringstream buf;

    buf << "N: " << N << "\n\n";
    buf << "Detector_type: " << detectorType << "\n\n";
    buf << "Flux_type: " << fluxType << "\n";
    buf << "Flux_dir: " << (verticalFlux ? "vertical" : "isotropic") << "\n";

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
    } else {
        buf << "\n";
    }
    buf << "}\n\n";

    buf << "Particles: [";
    if (fluxType == "PLAW") {
        buf << "gamma";
    } else if (fluxType == "SEP") {
        buf << "proton";
    } else if (fluxType == "Galactic") {
        buf << ReadValue("particle:");
    } else {
        buf << "gamma, proton, electron, alpha";
    }
    buf << "]\n";

    buf << "Energies:\n{\n\t";
    if (fluxType == "PLAW" || fluxType == "COMP") {
        buf << "gamma: ";
    } else if (fluxType == "SEP") {
        buf << "proton: ";
    } else if (fluxType == "Galactic") {
        buf << ReadValue("particle:") << ": ";
    } else if (fluxType == "Uniform") {
        buf << "gamma: (0.001, 500),\n\t";
        buf << "electron: (0.001, 100),\n\t";
        buf << "proton: (1, 10000),\n\t";
        buf << "alpha: (10, 10000),\n";
    }
    if (fluxType != "Uniform") {
        buf << "(" << ReadValue("E_min:") << ", " << ReadValue("E_max:") << ")\n";
    }
    buf << "}\n\n";

    buf << "Geometry:\n{\n\t";
    buf << "TunaCan: " << geometry->sizes.tunaCanThick / mm << "\n\t";
    buf << "Shell: " << sizes.shellThick / mm << "\n\t";
    buf << "Tyvek: " << sizes.tyvekThick / mm << "\n\t";
    buf << "Veto: " << sizes.vetoThick / mm << "\n\t";
    buf << "Gap: " << sizes.gapSize / mm << "\n\t";
    buf << "LED: " << sizes.LEDSize / mm << "\n}\n\n";

    buf << "Model_size:\n{\n\t";
    buf << "R: " << geometry->modelSize.y() / mm << "\n\t";
    buf << "H: " << geometry->modelSize.z() * 2. / mm << "\n}\n\n";

    buf << "Counts:\n{\n\t";
    buf << "Crystal_only: " << crystalOnly << "\n\t";
    buf << "Veto_then_Crystal: " << crystalAndVeto << "\n}\n\n";

    buf << "Rates:\n{\n\t";
    buf << std::fixed << std::setprecision(6);
    if (rate_ok) {
        buf << "Area: " << A_eff_cm2 << "\n\t";
        buf << "Integral: " << rr.integral << "\n\t";
        buf << "Ndot: " << rr.Ndot << "\n\t";
        buf << "Rate_Crystal_only: " << rr.rateCrystal << "\n\t";
        buf << "Rate_Both: " << rr.rateBoth << "\n";
    } else {
        buf << "Area: NaN\n\t";
        buf << "Integral: NaN\n\t";
        buf << "Ndot: NaN\n\t";
        buf << "Rate_Crystal_only: NaN\n\t";
        buf << "Rate_Both: NaN\n";
    }
    buf << "}\n\n";

    auto sanitize = [](std::string ss) {
        for (char &c: ss) if (c == ' ') c = '_';
        return ss;
    };

    std::string filename = "info_" + detectorType + "_" + fluxType;
    if (fluxType == "Galactic") {
        const std::string part = ReadValue("particle:");
        const std::string phi = ReadValue("phiMV:");
        filename += "_particle:" + part + "_phiMV:" + phi + ".txt";
    // } else if (fluxType == "SEP") {
    //     const std::string y = ReadValue("year:");
    //     const std::string order = ReadValue("order:");
    //     filename += "_year:" + y + "_order:" + order + ".txt";
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
