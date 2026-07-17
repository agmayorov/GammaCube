#include "Loader.hh"

#include <fstream>
#include <iostream>
#include <algorithm>
#include <cctype>
#include <ctime>
#include <sstream>
#include <unistd.h>

#ifdef G4UI_USE_QT
#include "G4UIQt.hh"
#include <QTimer>
#include <QApplication>
#include <QMainWindow>
#include <QDockWidget>
#include <QTabWidget>
#include <QTabBar>
#include <QToolBar>
#include <QMenuBar>
#include <QEventLoop>
#include <QAction>
#endif

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
    nBins = 1;
    saveSecondaries = false;
    savePhotons = false;

    bool visOnly = false;
    std::string visOutputFile = "frame_output";
    int visExportDelayMs = 3000;
    int visWindowWidth = 2000;
    int visWindowHeight = 3000;
    bool visFullScreen = false;

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
            if (Sizes::vetoTopRoundedRadius > 0) {
                Sizes::vetoChamferHeight = 0 * mm;
            }
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
        } else if (std::string(argv[i]) == "--frame-csv") {
            LoadFromCSV(argv[i + 1]);
        } else if (std::string input(argv[i]); input == "-v" || input == "--vis-only") {
            useUI = false;
            visOnly = true;
            macroFile = "../vis.mac";
        } else if (input == "--vis-output") {
            visOutputFile = argv[i + 1];
        } else if (input == "--vis-delay-ms") {
            visExportDelayMs = std::stoi(argv[i + 1]);
        } else if (input == "--vis-width") {
            visWindowWidth = std::stoi(argv[i + 1]);
        } else if (input == "--vis-height") {
            visWindowHeight = std::stoi(argv[i + 1]);
        } else if (input == "--vis-fullscreen") {
            visFullScreen = true;
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

        physicsList->RegisterPhysics(opticalPhysics);
    }

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

    area = Area_cm2(Sizes::modelRadius, Sizes::modelHeight, dir);

    runManager->SetUserInitialization(new ActionInitialization(area));
    runManager->Initialize();

    visManager = new G4VisExecutive;
    visManager->Initialize();

    G4UImanager* UImanager = G4UImanager::GetUIpointer();

    if (visOnly) {
#ifdef G4UI_USE_QT
        auto* ui = new G4UIExecutive(argc, argv, "qt");

        auto hideGeant4QtPanels = [ui]() {
            auto* qtSession = dynamic_cast<G4UIQt*>(ui->GetSession());
            if (!qtSession) {
                return;
            }

            if (auto* coutDock = qtSession->GetCoutDockWidget()) {
                coutDock->hide();

                if (coutDock->toggleViewAction()) {
                    coutDock->toggleViewAction()->setVisible(false);
                }
            }

            if (auto* uiDock = qtSession->GetUserInterfaceWidget()) {
                uiDock->hide();
            }

            if (QMainWindow* mainWindow = qtSession->GetMainWindow()) {
                for (auto* dock : mainWindow->findChildren<QDockWidget*>()) {
                    dock->hide();

                    if (dock->toggleViewAction()) {
                        dock->toggleViewAction()->setVisible(false);
                    }
                }

                for (auto* toolbar : mainWindow->findChildren<QToolBar*>()) {
                    toolbar->hide();
                }

                if (mainWindow->menuBar()) {
                    mainWindow->menuBar()->hide();
                }
            }

            QApplication::processEvents(QEventLoop::AllEvents, 100);
        };

        hideGeant4QtPanels();

        QTimer::singleShot(0, hideGeant4QtPanels);
        QTimer::singleShot(100, hideGeant4QtPanels);
        QTimer::singleShot(300, hideGeant4QtPanels);

        auto trimLocal = [](std::string s) {
            auto notSpace = [](unsigned char c) {
                return !std::isspace(c);
            };

            s.erase(s.begin(), std::find_if(s.begin(), s.end(), notSpace));
            s.erase(std::find_if(s.rbegin(), s.rend(), notSpace).base(), s.end());

            return s;
        };

        auto startsWith = [](const std::string& s, const std::string& prefix) {
            return s.rfind(prefix, 0) == 0;
        };

        auto executeVisMacroWithoutPrematureExit = [&]() {
            hideGeant4QtPanels();

            std::ifstream macro(macroFile);
            if (!macro.is_open()) {
                G4Exception(
                    "Loader::Loader",
                    "VIS_MACRO_OPEN_FAIL",
                    FatalException,
                    ("Cannot open visualisation macro: " + macroFile).c_str()
                );
            }

            G4cout << "[vis-only] Executing visualisation macro without exit/export commands: "
                   << macroFile << G4endl;

            std::string line;
            while (std::getline(macro, line)) {
                std::string cmd = trimLocal(line);

                if (cmd.empty()) {
                    continue;
                }

                if (cmd[0] == '#') {
                    continue;
                }

                if (cmd == "exit" ||
                    cmd == "/exit" ||
                    cmd == "/control/exit" ||
                    cmd == "/session/terminate" ||
                    startsWith(cmd, "/vis/ogl/export") ||
                    startsWith(cmd, "/vis/ogl/set/exportFormat")) {
                    G4cout << "[vis-only] skipped: " << cmd << G4endl;
                    continue;
                }

                G4int status = UImanager->ApplyCommand(cmd);

                if (status != 0) {
                    G4cerr << "[vis-only] command failed with status "
                           << status << ": " << cmd << G4endl;
                }

                hideGeant4QtPanels();
            }

            hideGeant4QtPanels();
        };

        auto prepareQtViewerForExport = [
            ui,
            hideGeant4QtPanels,
            visWindowWidth,
            visWindowHeight,
            visFullScreen
        ]() {
            hideGeant4QtPanels();

            auto* qtSession = dynamic_cast<G4UIQt*>(ui->GetSession());
            if (!qtSession) {
                G4cerr << "[vis-only] ERROR: current UI session is not G4UIQt." << G4endl;
                return;
            }

            QMainWindow* mainWindow = qtSession->GetMainWindow();
            if (!mainWindow) {
                G4cerr << "[vis-only] ERROR: G4UIQt main window is null." << G4endl;
                return;
            }

            if (auto* viewerTabs = qtSession->GetViewerTabWidget()) {
                if (auto* tabBar = viewerTabs->findChild<QTabBar*>()) {
                    tabBar->hide();
                }

                viewerTabs->setMinimumSize(visWindowWidth, visWindowHeight);
                viewerTabs->resize(visWindowWidth, visWindowHeight);
            }

            mainWindow->setMinimumSize(visWindowWidth, visWindowHeight);
            mainWindow->resize(visWindowWidth, visWindowHeight);

            if (visFullScreen) {
                mainWindow->showFullScreen();
            } else {
                mainWindow->showNormal();
                mainWindow->resize(visWindowWidth, visWindowHeight);
            }

            hideGeant4QtPanels();

            mainWindow->raise();
            mainWindow->activateWindow();

            QApplication::processEvents(QEventLoop::AllEvents, 500);

            hideGeant4QtPanels();

            G4cout << "[vis-only] Qt viewer prepared. Requested window size: "
                   << visWindowWidth << "x" << visWindowHeight
                   << ", fullscreen: " << (visFullScreen ? "true" : "false") << G4endl;
        };

        char cwdBuffer[4096];

        if (getcwd(cwdBuffer, sizeof(cwdBuffer)) != nullptr) {
            G4cout << "[vis-only] Current working directory: " << cwdBuffer << G4endl;
        }

        G4cout << "[vis-only] Export target basename: " << visOutputFile << G4endl;
        G4cout << "[vis-only] Export delay: " << visExportDelayMs << " ms" << G4endl;
        G4cout << "[vis-only] Requested Qt window size: "
               << visWindowWidth << "x" << visWindowHeight << G4endl;
        G4cout << "[vis-only] Fullscreen: "
               << (visFullScreen ? "true" : "false") << G4endl;

        executeVisMacroWithoutPrematureExit();
        hideGeant4QtPanels();

        QTimer::singleShot(0, hideGeant4QtPanels);
        QTimer::singleShot(100, hideGeant4QtPanels);
        QTimer::singleShot(300, hideGeant4QtPanels);

        QTimer::singleShot(
            visExportDelayMs,
            [UImanager, visOutputFile, prepareQtViewerForExport, hideGeant4QtPanels]() {
                G4cout << "[vis-only] Preparing viewer before PNG export..." << G4endl;

                hideGeant4QtPanels();
                prepareQtViewerForExport();
                hideGeant4QtPanels();

                UImanager->ApplyCommand("/vis/viewer/rebuild");
                UImanager->ApplyCommand("/vis/viewer/refresh");
                UImanager->ApplyCommand("/vis/viewer/update");
                UImanager->ApplyCommand("/vis/viewer/flush");

                hideGeant4QtPanels();

                QTimer::singleShot(
                    1000,
                    [UImanager, visOutputFile, hideGeant4QtPanels]() {
                        hideGeant4QtPanels();

                        G4cout << "[vis-only] Exporting PNG..." << G4endl;

                        G4int fmtStatus = UImanager->ApplyCommand("/vis/ogl/set/exportFormat png");
                        G4cout << "[vis-only] set exportFormat status: "
                               << fmtStatus << G4endl;

                        G4int exportStatus = UImanager->ApplyCommand(
                            G4String("/vis/ogl/export ") + visOutputFile
                        );

                        G4cout << "[vis-only] export status: "
                               << exportStatus << G4endl;

                        UImanager->ApplyCommand(
                            G4String("/control/shell ls -lh ") + visOutputFile + "*"
                        );

                        hideGeant4QtPanels();

                        QTimer::singleShot(
                            1500,
                            []() {
                                G4cout << "[vis-only] Closing Qt session." << G4endl;

                                if (qApp) {
                                    qApp->quit();
                                }
                            }
                        );
                    }
                );
            }
        );

        ui->SessionStart();

        delete ui;
#else
        G4Exception(
            "Loader::Loader",
            "QT_NOT_AVAILABLE",
            FatalException,
            "--vis-only requires Geant4 built with Qt support."
        );
#endif
    } else if (!useUI) {
        const G4String command = "/control/execute ";
        UImanager->ApplyCommand(command + macroFile);
    } else {
        auto* ui = new G4UIExecutive(argc, argv, "qt");

#ifdef G4UI_USE_QT
        auto hideOutputDock = [ui]() {
            auto* qtSession = dynamic_cast<G4UIQt*>(ui->GetSession());
            if (!qtSession) {
                return;
            }

            if (auto* coutDock = qtSession->GetCoutDockWidget()) {
                coutDock->hide();

                if (coutDock->toggleViewAction()) {
                    coutDock->toggleViewAction()->setVisible(false);
                }
            }

            QApplication::processEvents(QEventLoop::AllEvents, 100);
        };

        hideOutputDock();

        QTimer::singleShot(0, hideOutputDock);
        QTimer::singleShot(100, hideOutputDock);
        QTimer::singleShot(300, hideOutputDock);
#endif

        UImanager->ApplyCommand("/control/execute ../vis.mac");

#ifdef G4UI_USE_QT
        hideOutputDock();
#endif

        ui->SessionStart();

        delete ui;
    }
}

Loader::~Loader() {
    delete runManager;
    delete visManager;
}

std::string Loader::ReadValue(const std::string& key, const std::string& filepath) const {
    std::ifstream file(filepath.empty() ? configPath : filepath);

    if (!file.is_open()) {
        G4Exception(
            "Loader::ReadValue",
            "FILE_OPEN_FAIL",
            FatalException,
            ("Cannot open " + (filepath.empty() ? configPath : filepath)).c_str()
        );
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

        if (!token.empty()) {
            result.push_back(token);
        }
    }

    return result;
}

void Loader::LoadFromCSV(const std::string& path) {
    std::ifstream file(path);

    if (!file.is_open()) {
        G4Exception(
            "Loader::LoadFromCSV",
            "FILE_OPEN_FAIL",
            FatalException,
            ("Cannot open " + path).c_str()
        );
    }

    std::string line;
    std::getline(file, line);

    while (std::getline(file, line)) {
        auto tokens = Split(line);

        if (tokens.size() < 3) {
            continue;
        }

        std::string name = tokens[0];
        double deg = std::stod(tokens[1]);
        double shift = std::stod(tokens[2]);

        if (name == "tyvekOut") {
            DEG::tyvekOut = deg;
            shift::tyvekOut = shift;
        } else if (name == "veto") {
            DEG::veto = deg;
            shift::veto = shift;
        } else if (name == "tyvekMid") {
            DEG::tyvekMid = deg;
            shift::tyvekMid = shift;
        } else if (name == "opticLayerVeto") {
            DEG::opticLayerVeto = deg;
            shift::opticLayerVeto = shift;
        } else if (name == "rubber") {
            DEG::rubber = deg;
            shift::rubber = shift;
        } else if (name == "shell") {
            DEG::shell = deg;
            shift::shell = shift;
        } else if (name == "bottomVetoShell") {
            DEG::bottomVetoShell = deg;
            shift::bottomVetoShell = shift;
        } else if (name == "tyvekBottom") {
            DEG::tyvekBottom = deg;
            shift::tyvekBottom = shift;
        } else if (name == "bottomVeto") {
            DEG::bottomVeto = deg;
            shift::bottomVeto = shift;
        } else if (name == "opticLayerBottomVeto") {
            DEG::opticLayerBottomVeto = deg;
            shift::opticLayerBottomVeto = shift;
        } else if (name == "crystalContainer") {
            DEG::crystalContainer = deg;
            shift::crystalContainer = shift;
        } else if (name == "crystal") {
            DEG::crystal = deg;
            shift::crystal = shift;
        } else if (name == "tyvekIn") {
            DEG::tyvekIn = deg;
            shift::tyvekIn = shift;
        } else if (name == "crystalShell") {
            DEG::crystalShell = deg;
            shift::crystalShell = shift;
        } else if (name == "crystallGlass") {
            DEG::crystallGlass = deg;
            shift::crystallGlass = shift;
        } else if (name == "opticLayerCrystall") {
            DEG::opticLayerCrystall = deg;
            shift::opticLayerCrystall = shift;
        } else if (name == "holder") {
            DEG::holder = deg;
            shift::holder = shift;
        } else if (name == "crystalSiPM") {
            DEG::crystalSiPM = deg;
            shift::crystalSiPM = shift;
        } else if (name == "vetoSiPM") {
            DEG::vetoSiPM = deg;
            shift::vetoSiPM = shift;
        } else if (name == "vetoSpring") {
            DEG::vetoSpring = deg;
            shift::vetoSpring = shift;
        } else if (name == "vetoBoard") {
            DEG::vetoBoard = deg;
            shift::vetoBoard = shift;
        } else if (name == "bottomVetoSiPM") {
            DEG::bottomVetoSiPM = deg;
            shift::bottomVetoSiPM = shift;
        } else if (name == "tunaCan") {
            DEG::tunaCan = deg;
            shift::tunaCan = shift;
        } else if (name == "plate") {
            DEG::plate = deg;
            shift::plate = shift;
        } else if (name == "plateHole") {
            DEG::plateHole = deg;
            shift::plateHole = shift;
        }
    }

    file.close();
}