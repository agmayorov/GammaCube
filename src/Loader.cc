#include "Loader.hh"


Loader::Loader(int argc, char **argv) {
    numThreads = G4Threading::G4GetNumberOfCores();
    useUI = true;
    macroFile = "../run.mac";
    temperature = 20;
    detectorType = "NaI";
    sizes.vetoThick = 0 * mm;
    sizes.shellThick = 0 * mm;
    sizes.gapSize = 0 * mm;
    sizes.tapeThick = 1 * mm;
    sizes.LEDSize = 5 * mm;
    verticalFlux = false;
    doubleLED = false;

    for (int i = 0; i < argc; i++) {
        std::string input(argv[i]);
        if (input == "-i" || input == "-input") {
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
        } else if (input == "-tape") {
            sizes.tapeThick = std::stod(argv[i + 1]) * mm;
        } else if (input == "--double-LED") {
            doubleLED = true;
        } else if (input == "-LED") {
            sizes.LEDSize = std::stod(argv[i + 1]);
        } else if (input == "-d" || input == "--detector") {
            detectorType = argv[i + 1];
        } else if (input == "--vertical-flux") {
            verticalFlux = true;
        }
    }

    CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);
    CLHEP::HepRandom::setTheSeed(time(NULL));

#ifdef G4MULTITHREADED
    runManager = new G4MTRunManager;
    runManager->SetNumberOfThreads(runManager->GetNumberOfThreads());
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

    runManager->SetUserInitialization(new ActionInitialization(realWorld, verticalFlux));

    visManager = new G4VisExecutive;
    visManager->Initialize();
    G4UImanager *UImanager = G4UImanager::GetUIpointer();

    if (!useUI) {
        G4String command = "/control/execute ";
        UImanager->ApplyCommand(command + macroFile);
    } else {
        G4UIExecutive *ui = new G4UIExecutive(argc, argv, "qt");
        UImanager->ApplyCommand("/control/execute ../vis.mac");
        ui->SessionStart();
        delete ui;
    }
}

Loader::~Loader() {
    delete runManager;
    delete visManager;
}
