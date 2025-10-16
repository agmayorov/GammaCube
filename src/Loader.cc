#include "Loader.hh"

Loader::Loader(int argc, char **argv) {
    numThreads = G4Threading::G4GetNumberOfCores();
    useUI = true;
    macroFile = "../run.mac";
    detectorType = "CsI";
    fluxType = "Uniform";
    geomConfigPath = "../geometry_config.txt";
    verticalFlux = false;

    for (int i = 0; i < argc; i++) {
        if (std::string input(argv[i]); input == "-i" || input == "-input") {
            macroFile = argv[i + 1];
            useUI = false;
        } else if (input == "-t" || input == "-threads") {
            numThreads = std::stoi(argv[i + 1]);
        } else if (input == "-noUI") {
            useUI = false;
        } else if (input == "-d" || input == "--detector") {
            detectorType = argv[i + 1];
        } else if (input == "-f" || input == "--flux-type") {
            fluxType = argv[i + 1];
        } else if (input == "--vertical-flux") {
            verticalFlux = true;
        } else if (input == "-g" || input == "--geom-config") {
            geomConfigPath = argv[i + 1];
        }
    }

    configPath = "../Flux_config/" + fluxType + "_params.txt";
    ParseGeomConfig();

    CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);
    CLHEP::HepRandom::setTheSeed(time(nullptr));

#ifdef G4MULTITHREADED
    runManager = new G4MTRunManager;
    runManager->SetNumberOfThreads(numThreads);
#else
    runManager = new G4RunManager;
#endif

    Geometry *realWorld = new Geometry(detectorType);
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

    runManager->SetUserInitialization(new ActionInitialization(verticalFlux, fluxType));
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
    SaveConfig();
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


void Loader::SaveConfig() const {
    const int N = std::stoi(ReadValue("/run/beamOn", "../run.mac"));

    const FluxDir dir = verticalFlux ? FluxDir::Vertical : FluxDir::Isotropic;
    double A_eff_cm2 = effectiveArea_cm2(Sizes::modelRadius, Sizes::modelHeight, dir);

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


inline std::string trim(std::string st) {
    auto notspace = [](const unsigned char c) { return !std::isspace(c); };
    st.erase(st.begin(), std::find_if(st.begin(), st.end(), notspace));
    st.erase(std::find_if(st.rbegin(), st.rend(), notspace).base(), st.end());
    return st;
}

inline std::string to_upper(std::string st) {
    std::transform(st.begin(), st.end(), st.begin(),
                   [](const unsigned char c) { return std::toupper(c); });
    return st;
}

void Loader::ParseGeomConfig() {
    std::regex kvRe(R"(^\s*([A-Za-z0-9_]+)\s*:\s*(.*?)\s*$)");

    auto setD = [](G4double &ref, const char *keyName) {
        return [&ref, keyName](const std::string &v) {
            try { ref = std::stod(v) * mm; } catch (...) {
                std::cerr << "[WARN] Bad double: " << keyName << ": \"" << v << "\"\n";
            }
        };
    };

    auto setI = [](G4int &ref, const char *keyName) {
        return [&ref, keyName](const std::string &v) {
            try { ref = std::stoi(v); } catch (...) {
                std::cerr << "[WARN] Bad int: " << keyName << ": \"" << v << "\"\n";
            }
        };
    };

    auto stripInlineComment = [](std::string &str) {
        auto pos = str.find('#');
        if (pos != std::string::npos) str.erase(pos);
        str = trim(str);
    };

    std::unordered_map<std::string, std::function<void(const std::string &)> > setters = {
        {"Model_height", setD(Sizes::modelHeight, "Model_height")},
        {"Model_radius", setD(Sizes::modelRadius, "Model_radius")},

        {"TunaCan_thick_wall", setD(Sizes::tunaCanThickWall, "TunaCan_thick_wall")},
        {"TunaCan_thick_top", setD(Sizes::tunaCanThickTop, "TunaCan_thick_top")},
        {"TunaCan_thick_bottom", setD(Sizes::tunaCanThickBottom, "TunaCan_thick_bottom")},

        {"Gap_size_wall", setD(Sizes::gapSizeWall, "Gap_size_wall")},
        {"Gap_size_top", setD(Sizes::gapSizeTop, "Gap_size_top")},
        {"Gap_size_bottom", setD(Sizes::gapSizeBottom, "Gap_size_bottom")},

        {"TyvekOut_thick_wall", setD(Sizes::tyvekOutThickWall, "TyvekOut_thick_wall")},
        {"TyvekOut_thick_top", setD(Sizes::tyvekOutThickTop, "TyvekOut_thick_top")},
        {"TyvekOut_thick_bottom", setD(Sizes::tyvekOutThickBottom, "TyvekOut_thick_bottom")},

        {"Veto_thick_wall", setD(Sizes::vetoThickWall, "Veto_thick_wall")},
        {"Veto_thick_top", setD(Sizes::vetoThickTop, "Veto_thick_top")},
        {"Veto_thick_bottom", setD(Sizes::vetoThickBottom, "Veto_thick_bottom")},

        {"Veto_chamfer_heigh", setD(Sizes::vetoChamferHeigh, "Veto_chamfer_heigh")},

        {"TyvekMid_thick_wall", setD(Sizes::tyvekMidThickWall, "TyvekMid_thick_wall")},
        {"TyvekMid_thick_top", setD(Sizes::tyvekMidThickTop, "TyvekMid_thick_top")},
        {"TyvekMid_thick_bottom", setD(Sizes::tyvekMidThickBottom, "TyvekMid_thick_bottom")},

        {"Rubber_radius", setD(Sizes::rubberRadius, "Rubber_radius")},
        {"Rubber_height", setD(Sizes::rubberHeight, "Rubber_height")},

        {"Aluminium_thick_wall", setD(Sizes::AlThickWall, "Aluminium_thick_wall")},
        {"Aluminium_thick_top", setD(Sizes::AlThickTop, "Aluminium_thick_top")},

        {"AluminiumCap_thick_wall", setD(Sizes::AlCapThickWall, "AluminiumCap_thick_wall")},
        {"AluminiumCap_thick_bottom", setD(Sizes::AlCapThickBottom, "AluminiumCap_thick_bottom")},

        {"TyvekIn_thick_wall", setD(Sizes::tyvekInThickWall, "TyvekIn_thick_wall")},
        {"TyvekIn_thick_top", setD(Sizes::tyvekInThickTop, "TyvekIn_thick_top")},
        {"TyvekIn_thick_bottom", setD(Sizes::tyvekInThickBottom, "TyvekIn_thick_bottom")},

        {"CrystalLED_count", setI(Sizes::crystalLEDCount, "CrystalLED_count")},
        {"CrystalLED_radius", setD(Sizes::crystalLEDRadius, "CrystalLED_radius")},
        {"CrystalLED_height", setD(Sizes::crystalLEDHeight, "CrystalLED_height")},

        {"VetoLED_count", setI(Sizes::vetoLEDCount, "VetoLED_count")},
        {"VetoLED_radius", setD(Sizes::vetoLEDRadius, "VetoLED_radius")},
        {"VetoLED_height", setD(Sizes::vetoLEDHeight, "VetoLED_height")},

        {"VetoLED_bottom_count", setI(Sizes::vetoBottomLEDCount, "VetoLED_bottom_count")},
        {"VetoLED_bottom_radius", setD(Sizes::vetoBottomLEDRadius, "VetoLED_bottom_radius")},
        {"VetoLED_bottom_height", setD(Sizes::vetoBottomLEDHeight, "VetoLED_bottom_height")},

        {"Wire_radius", setD(Sizes::wireRadius, "Wire_radius")},
        {"Wire_insulation_thick", setD(Sizes::wireInsulationThick, "Wire_insulation_thick")},

        {"Pin_count", setI(Sizes::pinCount, "Pin_count")},
        {"Pin_radius", setD(Sizes::pinRadius, "Pin_radius")},

        {"Board_height", setD(Sizes::boardHeight, "Board_height")},
        {"Board_length", setD(Sizes::boardLength, "Board_length")},
        {"Board_width", setD(Sizes::boardWidth, "Board_width")},

        {"Board_space", setD(Sizes::boardSpace, "Board_space")},
    };

    std::string line;
    size_t lineNo = 0;

    std::ifstream geomConfigFile(geomConfigPath);
    while (std::getline(geomConfigFile, line)) {
        ++lineNo;

        std::string trimmed = trim(line);
        if (trimmed.empty() || (!trimmed.empty() && trimmed[0] == '#'))
            continue;

        std::smatch match;
        if (!std::regex_match(trimmed, match, kvRe)) {
            std::cerr << "[WARN] Line " << lineNo << ": cannot parse: " << line << "\n";
            continue;
        }

        std::string key = match[1].str();
        std::string val = trim(match[2].str());
        stripInlineComment(val);

        if (key == "Detector_type") {
            detectorType = val;
            continue;
        }

        auto it = setters.find(key);
        if (it == setters.end()) {
            std::cerr << "[WARN] Line " << lineNo << ": unknown key \"" << key << "\"\n";
            continue;
        }

        it->second(val);
    }

    Sizes::tunaCanAllSize = Sizes::tunaCanThickTop + Sizes::tunaCanThickBottom + Sizes::tunaCanThickWall;
    Sizes::tyvekInAllSize = Sizes::tyvekInThickTop + Sizes::tyvekInThickBottom + Sizes::tyvekInThickWall;
    Sizes::tyvekOutAllSize = Sizes::tyvekOutThickTop + Sizes::tyvekOutThickBottom + Sizes::tyvekOutThickWall;
    Sizes::tyvekMidAllSize = Sizes::tyvekMidThickTop + Sizes::tyvekMidThickBottom + Sizes::tyvekMidThickWall;
    Sizes::vetoAllSize = Sizes::vetoThickTop + Sizes::vetoThickBottom + Sizes::vetoThickWall;
    Sizes::AlAllSize = Sizes::AlThickTop + Sizes::AlCapThickBottom + Sizes::AlThickWall + Sizes::AlCapThickWall;
    Sizes::rubberAllSize = Sizes::rubberRadius + Sizes::rubberHeight;
    Sizes::crystalLEDAllSize = Sizes::crystalLEDCount + Sizes::crystalLEDRadius + Sizes::crystalLEDHeight;
    Sizes::vetoLEDAllSize = Sizes::vetoLEDCount + Sizes::vetoLEDRadius + Sizes::vetoLEDHeight;
    Sizes::vetoBottomLEDAllSize = Sizes::vetoBottomLEDCount + Sizes::vetoBottomLEDRadius + Sizes::vetoBottomLEDHeight;
}
