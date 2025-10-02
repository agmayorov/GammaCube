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
    const std::string filepath = "info.txt";
    std::ofstream out(filepath);
    if (!out.is_open()) {
        G4cerr << "Ошибка: не удалось открыть файл " << filepath << G4endl;
        return;
    }

    out << "N: " << std::stoi(ReadValue("/run/beamOn", "../run.mac")) << "\n\n";
    out << "Detector_type: " << detectorType << "\n\n";
    out << "Flux_type: " << fluxType << "\n";
    out << "Flux_dir: " << (verticalFlux ? "vertical" : "isotropic") << "\n";

    out << "Flux_params:\n{\n\t";
    if (fluxType == "PLAW") {
        out << "A: " << std::stod(ReadValue("A:")) << ",\n\t";
        out << "alpha: " << std::stod(ReadValue("alpha:")) << ",\n\t";
        out << "E_Piv: " << std::stod(ReadValue("E_Piv:")) << "\n";
    } else if (fluxType == "COMP") {
        out << "A: " << std::stod(ReadValue("A:")) << ",\n\t";
        out << "alpha: " << std::stod(ReadValue("alpha:")) << ",\n\t";
        out << "E_Piv: " << std::stod(ReadValue("E_Piv:")) << ",\n\t";
        out << "E_Peak: " << std::stod(ReadValue("E_Piv:")) << "\n";
    } else if (fluxType == "SEP") {
        out << "year: " << std::stoi(ReadValue("year:")) << ",\n\t";
        out << "order: " << std::stoi(ReadValue("order:")) << "\n";
    } else if (fluxType == "Galactic") {
        out << "phiMV: " << std::stoi(ReadValue("phiMV:")) << "\n";
    } else {
        out << "\n";
    }
    out << "}\n\n";

    out << "Particles: [";
    if (fluxType == "PLAW") {
        out << "gamma";
    } else if (fluxType == "SEP" or fluxType == "Galactic") {
        out << "proton";
    } else {
        out << "gamma, proton, electron, alpha";
    }
    out << "]\n";

    out << "Energies:\n{\n\t";
    if (fluxType == "PLAW" or fluxType == "COMP") {
        out << "gamma: ";
    } else if (fluxType == "SEP") {
        out << "proton: ";
    } else if (fluxType == "Galactic") {
        out << ReadValue("particle:") << ": ";
    } else if (fluxType == "Uniform") {
        out << "gamma: (0.001, 500),\n\t";
        out << "electron: (0.001, 100),\n\t";
        out << "proton: (1, 10000),\n\t";
        out << "alpha: (10, 10000),\n";
    }
    if (fluxType != "Uniform") {
        out << "(" << ReadValue("E_min:") << ", ";
        out << ReadValue("E_max:") << ")\n";
    }
    out << "}\n\n";

    out << "Geometry:\n{\n\t";
    out << "TunaCan: " << geometry->sizes.tunaCanThick / mm << "\n\t";
    out << "Shell: " << sizes.shellThick / mm << "\n\t";
    out << "Tyvek: " << sizes.tyvekThick / mm << "\n\t";
    out << "Veto: " << sizes.vetoThick / mm << "\n\t";
    out << "Gap: " << sizes.gapSize / mm << "\n\t";
    out << "LED: " << sizes.LEDSize / mm << "\n}\n\n";

    out << "Model_size:\n{\n\t";
    out << "R: " << geometry->modelSize.y() / mm << "\n\t";
    out << "H: " << geometry->modelSize.z() * 2. / mm << "\n}";

    out.close();
    std::cout << "Configuration saved in " << filepath << std::endl;
}
