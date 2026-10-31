#include "Loader.hh"

#include <cctype>

using namespace Configuration;

namespace {
    void ArgumentError(const std::string& message) {
        G4Exception("Loader::Loader", "Arguments", FatalException, message.c_str());
    }

    // A flag, not a value: negative numbers such as "--source-z -10" are values
    bool IsFlag(const std::string& token) {
        return token.size() > 1 && token[0] == '-' && !std::isdigit(static_cast<unsigned char>(token[1])) &&
               token[1] != '.';
    }

    // The whole value must be a number of this type: "12abc" and, for an integer, "8.5" are errors
    template <typename T>
    T Number(const std::string& flag, const std::string& text) {
        std::istringstream in(text);
        T number{};
        if (!(in >> number) || !in.eof()) ArgumentError("Flag " + flag + " has an invalid value: " + text);
        return number;
    }

    // -c/--config file: the same flags as on the command line, split on whitespace; '#' starts a comment
    void ReadConfigFile(const std::string& path, std::vector<std::string>& out) {
        std::ifstream file(path);
        if (!file.is_open()) ArgumentError("Cannot open config file: " + path);

        std::string line;
        while (std::getline(file, line)) {
            std::istringstream tokens(line.substr(0, line.find('#')));
            for (std::string token; tokens >> token;) {
                if (token == "-c" || token == "--config") ArgumentError("Config files cannot include others: " + path);
                out.push_back(token);
            }
        }
    }

    // Config files go first, so that any flag given on the command line overrides them wherever -c stands
    std::vector<std::string> CollectArguments(const int argc, char** argv) {
        std::vector<std::string> fromFiles;
        std::vector<std::string> fromCommandLine;
        for (int i = 1; i < argc; ++i) {
            const std::string token(argv[i]);
            if (token != "-c" && token != "--config") {
                fromCommandLine.push_back(token);
                continue;
            }
            if (i + 1 >= argc || IsFlag(argv[i + 1])) ArgumentError("Flag " + token + " needs a value");
            ReadConfigFile(argv[++i], fromFiles);
        }
        fromFiles.insert(fromFiles.end(), fromCommandLine.begin(), fromCommandLine.end());
        return fromFiles;
    }
}

Loader::Loader(int argc, char** argv) {
    numThreads = G4Threading::G4GetNumberOfCores();
    useUI = true;
    macroFile = "../run.mac";
    detectorType = "CsI";
    fluxType = "Uniform";
    fluxDirection = "isotropic";
    eCrystalThreshold = 0 * MeV;
    eVetoThreshold = 0 * MeV;
    useOptics = false;
    yieldScale = 1;
    oCrystalThreshold = 0 * MeV;
    oVetoThreshold = 0 * MeV;
    photonTimeLimit = -1.;
    sipmCells = 0;
    sipmRecoveryTime = 95 * ns;
    sipmCrosstalk = 0.;
    outputFile = "GammaCube.root";
    nBins = 1000;
    saveSecondaries = false;
    savePhotons = false;

    const std::vector<std::string> args = CollectArguments(argc, argv);
    for (const auto& arg : args) runArguments += (runArguments.empty() ? "" : " ") + arg;

    for (size_t i = 0; i < args.size(); i++) {
        const std::string& input = args[i];
        auto value = [&]() -> std::string {
            if (i + 1 >= args.size() || IsFlag(args[i + 1])) ArgumentError("Flag " + input + " needs a value");
            return args[++i];
        };
        auto real = [&]() { return Number<G4double>(input, value()); };
        auto integer = [&]() { return Number<G4int>(input, value()); };

        if (input == "-i" || input == "--input") {
            macroFile = value();
            useUI = false;
            viewDeg = 360 * deg;
        } else if (input == "-t" || input == "--threads") {
            numThreads = integer();
        } else if (input == "-ys" || input == "--yield-scale") {
            yieldScale = integer();
        } else if (input == "--bins") {
            nBins = integer();
        } else if (input == "-vd" || input == "--view-deg") {
            const G4double angle = real() * deg;
            if (useUI) viewDeg = angle;
        } else if (input == "-ct" || input == "--crystal-threshold") {
            eCrystalThreshold = real() * MeV;
        } else if (input == "-vt" || input == "--veto-threshold") {
            eVetoThreshold = real() * MeV;
        } else if (input == "-oct" || input == "--crystal-optic-threshold") {
            oCrystalThreshold = integer();
        } else if (input == "-ovt" || input == "--veto-optic-threshold") {
            oVetoThreshold = integer();
        } else if (input == "-obvt" || input == "--bottom-veto-optic-threshold") {
            oBottomVetoThreshold = integer();
        } else if (input == "-ptl" || input == "--photon-time-limit") {
            photonTimeLimit = real() * us;
        } else if (input == "-sc" || input == "--sipm-cells") {
            sipmCells = integer();
        } else if (input == "-sr" || input == "--sipm-recovery") {
            sipmRecoveryTime = real() * ns;
        } else if (input == "-xt" || input == "--sipm-crosstalk") {
            sipmCrosstalk = real();
        } else if (input == "-vch" || input == "--veto-chamfer-height") {
            Sizes::vetoChamferHeight = real() * mm;
        } else if (input == "-vtr" || input == "--veto-top-rounded") {
            Sizes::vetoTopRoundedRadius = real() * mm;
        } else if (input == "-noUI") {
            useUI = false;
        } else if (input == "--polished") {
            polishedTyvek = true;
        } else if (input == "-d" || input == "--detector") {
            detectorType = value();
        } else if (input == "-csc" || input == "--crystal-sipm-config" || input == "-sipm") {
            crystalSiPMConfig = value();
        } else if (input == "-f" || input == "--flux-type") {
            fluxType = value();
        } else if (input == "--flux-dir" || input == "--f-dir" || input == "-fd") {
            fluxDirection = value();
        } else if (input == "--use-optics") {
            useOptics = true;
        } else if (input == "--save-secondaries") {
            saveSecondaries = true;
        } else if (input == "--save-photons") {
            savePhotons = true;
        } else if (input == "-gc" || input == "--geometry-config") {
            geometryConfig = value();
        } else if (input == "--engine") {
            engineConfig = value();
        } else if (input == "--camera") {
            // the model name is optional: "--camera" alone means the default camera
            const G4bool hasValue = i + 1 < args.size() && !IsFlag(args[i + 1]);
            cameraConfig = hasValue ? value() : "galper";
        } else if (input == "--trigger") {
            triggerType = value();
        } else if (input == "-w4" || input == "--weight4") {
            weight4 = real();
        } else if (input == "-w8" || input == "--weight8") {
            weight8 = real();
        } else if (input == "--source-r") {
            sourceR = real() * mm;
        } else if (input == "--source-z") {
            sourceZ = real() * mm;
        } else if (input == "--cone") {
            coneAngle = real() * deg;
        } else if (input == "--theta-min") {
            sphereThetaMin = real() * deg;
        } else if (input == "--theta-max") {
            sphereThetaMax = real() * deg;
        } else if (input == "--phi-min") {
            spherePhiMin = real() * deg;
        } else if (input == "--phi-max") {
            spherePhiMax = real() * deg;
        } else if (input == "--axis-theta") {
            sphereAxisTheta = real() * deg;
        } else if (input == "--axis-phi") {
            sphereAxisPhi = real() * deg;
        } else if (input == "--theta") {
            beamTheta = real() * deg;
        } else if (input == "--phi") {
            beamPhi = real() * deg;
        } else if (input == "-s" || input == "--seed") {
            seed = Number<G4long>(input, value());
        } else if (input == "-o" || input == "--output-file") {
            outputFile = value();
            outputFile += ".root";
        } else {
            ArgumentError("Unknown argument: " + input);
        }
    }

    // the rounded top replaces the chamfer
    if (Sizes::vetoTopRoundedRadius > 0) Sizes::vetoChamferHeight = 0 * mm;

    savePhotons = savePhotons and useOptics;
    if (!useOptics) {
        photonTimeLimit = -1.;
        sipmCells = 0;
        sipmCrosstalk = 0.;
    }

    if (fluxDirection == "vertical_down") {
        fluxDirection = "flat";
        beamTheta = 0 * deg;
        beamPhi = 0 * deg;
    } else if (fluxDirection == "vertical_up") {
        fluxDirection = "flat";
        beamTheta = 180 * deg;
        beamPhi = 0 * deg;
    } else if (fluxDirection == "horizontal") {
        fluxDirection = "flat";
        beamTheta = 90 * deg;
        beamPhi = 0 * deg;
    }
    if (triggerType != "gamma" && triggerType != "charged") {
        G4Exception("Loader::Loader", "Trigger", FatalException,
                    ("Trigger not found: " + triggerType + ".\nAvailable triggers: gamma, charged").c_str());
    }
    Layout::Configure();
    genSurface = GenSurface::For(fluxDirection);

    configPath = "../Flux_config/" + fluxType + "_params.txt";

    CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);
    if (seed == 0) seed = time(nullptr);
    CLHEP::HepRandom::setTheSeed(seed);

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

    if (fluxType == "Uniform") {
        std::string isLogStr = ReadValue("is_log:", configPath);
        isLogBin = isLogStr == "1" || isLogStr == "true";
    }

    area = genSurface.Norm_cm2();
    runManager->SetUserInitialization(new ActionInitialization(area));
    runManager->Initialize();

    visManager = new G4VisExecutive;
    visManager->Initialize();
    genSurfaceVis = new GenSurfaceVis;
    visManager->RegisterRunDurationUserVisAction("GenSurface", genSurfaceVis, GenSurfaceVis::Extent());
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
        killedPhotons = runAction->GetKilledPhotons();
        trackedPhotons = runAction->GetTrackedPhotons();
    }
    SaveConfig();
    RunPostProcessing();
}

Loader::~Loader() {
    delete runManager;
    delete visManager;
    delete genSurfaceVis;
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

    buf << "N: " << N << "\n";
    buf << "Seed: " << seed << "\n";
    buf << "Arguments: " << runArguments << "\n\n";
    buf << "Geometry_config: " << geometryConfig << "\n";
    buf << "Engine: " << (engineConfig.empty() ? "none" : engineConfig) << "\n";
    buf << "Camera: " << (cameraConfig.empty() ? "none" : cameraConfig) << "\n";
    buf << "Trigger: " << triggerType << (triggerType == "charged" ? " (Crystal && Veto)" : " (Crystal && !Veto)")
        << "\n\n";
    buf << "Detector_type: " << detectorType << "\n";
    buf << "Crystal_SiPM_configuration: " << crystalSiPMConfig << "\n";
    buf << "Tyvek_surface: " << (polishedTyvek ? "polished" : "diffuse") << "\n\n";
    buf << "Use_optics: " << useOptics << "\n";
    if (useOptics) {
        buf << "Photon_time_limit: ";
        if (photonTimeLimit > 0) {
            buf << photonTimeLimit / us << " us\n";
            buf << "Killed_photons: " << static_cast<long long>(killedPhotons)
                << " of " << static_cast<long long>(trackedPhotons);
            if (trackedPhotons > 0)
                buf << " (" << 100. * killedPhotons / trackedPhotons << " %)";
            buf << "\n";
        } else {
            buf << "off\n";
        }
        buf << "SiPM_microcells: ";
        if (sipmCells > 0)
            buf << sipmCells << " (recovery " << sipmRecoveryTime / ns << " ns)\n";
        else
            buf << "off\n";
        buf << "SiPM_crosstalk: ";
        if (sipmCrosstalk > 0) buf << sipmCrosstalk << "\n";
        else buf << "off\n";
    }
    buf << "\n";
    buf << "Flux_type: " << fluxType << "\n";
    buf << "Flux_dir: " << fluxDirection << "\n\n";

    buf << "Generation_surface:\n{\n\t";
    buf << "shape: " << genSurface.ShapeName() << ",\n\t";
    buf << "centre: (" << genSurface.Origin().x() / mm << ", " << genSurface.Origin().y() / mm << ", "
        << genSurface.Origin().z() / mm << ") mm,\n\t";
    if (genSurface.IsBox()) {
        buf << "box_half_size: (" << genSurface.BoxHalfSize().x() / mm << ", " << genSurface.BoxHalfSize().y() / mm
            << ", " << genSurface.BoxHalfSize().z() / mm << ") mm,\n\t";
    } else {
        buf << "cylinder_radius: " << genSurface.CylinderRadius() / mm << " mm,\n\t";
        buf << "cylinder_half_height: " << genSurface.CylinderHalfHeight() / mm << " mm,\n\t";
    }
    if (genSurface.IsIsotropic()) {
        buf << "R_gen: " << genSurface.Radius() / mm << " mm,\n\t";
        if (genSurface.IsSpherePart()) {
            buf << "theta_range: (" << genSurface.ThetaMin() / deg << ", " << genSurface.ThetaMax() / deg
                << ") deg,\n\t";
            buf << "phi_range: (" << genSurface.PhiMin() / deg << ", "
                << (genSurface.PhiMin() + genSurface.PhiWidth()) / deg << ") deg,\n\t";
            buf << "sphere_axis: (" << genSurface.SphereAxis().x() << ", " << genSurface.SphereAxis().y() << ", "
                << genSurface.SphereAxis().z() << "),\n\t";
            buf << "sphere_fraction: " << genSurface.SphereFraction() << ",\n\t";
        }
        buf << "geometric_factor: " << genSurface.GeomFactor_cm2sr() << " cm^2*sr\n";
    } else if (genSurface.IsPoint()) {
        buf << "theta: " << genSurface.Theta() / deg << " deg,\n\t";
        buf << "phi: " << genSurface.Phi() / deg << " deg,\n\t";
        buf << "axis: (" << genSurface.Axis().x() << ", " << genSurface.Axis().y() << ", "
            << genSurface.Axis().z() << "),\n\t";
        buf << "source_radius: " << genSurface.SourceRadius() / mm << " mm,\n\t";
        buf << "source_position: (" << genSurface.SourcePosition().x() / mm << ", "
            << genSurface.SourcePosition().y() / mm << ", " << genSurface.SourcePosition().z() / mm << ") mm";
        if (genSurface.IsCone()) {
            buf << ",\n\tcone_half_angle: " << genSurface.ConeAngle() / deg << " deg,\n\t";
            buf << "solid_angle_fraction: " << genSurface.ConeSolidAngleFraction();
        }
        buf << "\n";
    } else {
        buf << "theta: " << genSurface.Theta() / deg << " deg,\n\t";
        buf << "phi: " << genSurface.Phi() / deg << " deg,\n\t";
        buf << "axis: (" << genSurface.Axis().x() << ", " << genSurface.Axis().y() << ", "
            << genSurface.Axis().z() << "),\n\t";
        buf << "half_u: " << genSurface.HalfU() / mm << " mm,\n\t";
        buf << "half_v: " << genSurface.HalfV() / mm << " mm,\n\t";
        buf << "standoff: " << genSurface.Standoff() / mm << " mm,\n\t";
        buf << "S_perp: " << genSurface.SPerp_cm2() << " cm^2\n";
    }
    buf << "}\n\n";

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

    auto efficiencyBlock = [&](const std::string& title, const G4int only, const G4int withVeto) {
        const G4double eff = N > 0 ? static_cast<G4double>(only) / N : 0.;
        const G4double effBoth = N > 0 ? static_cast<G4double>(only + withVeto) / N : 0.;
        buf << title << ":\n{\n\t";
        buf << "Crystal_only: " << eff << "\n\t";
        buf << "Both: " << effBoth;
        if (genSurface.IsCone()) {
            buf << "\n\tSolid_angle_fraction: " << genSurface.ConeSolidAngleFraction() << "\n\t";
            buf << "Absolute_Crystal_only: " << eff * genSurface.ConeSolidAngleFraction() << "\n\t";
            buf << "Absolute_Both: " << effBoth * genSurface.ConeSolidAngleFraction();
        }
        buf << "\n}\n\n";
    };

    if (genSurface.IsPoint()) {
        buf << std::defaultfloat << std::setprecision(6);
        efficiencyBlock("Efficiency", crystalOnly, crystalAndVeto);
        if (useOptics) efficiencyBlock("Optical_efficiency", crystalOnlyOpt, crystalAndVetoOpt);
    } else {
        std::string area_dim = genSurface.IsIsotropic() ? " sr * cm^2" : " cm^2";
        std::string area_dim_inv = genSurface.IsIsotropic() ? " sr^-1 * cm^-2" : " cm^-2";

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
    }

    auto sanitize = [](std::string ss) {
        for (char& c : ss) if (c == ' ') c = '_';
        return ss;
    };

    std::string filename = "info_" + detectorType + "_" + fluxType;
    if (fluxType == "Galactic") {
        const std::string part = ReadValue("particle:");
        const std::string phi = ReadValue("phiMV:");
        filename += "_particle:" + part + "_phiMV:" + phi;
    } else if (fluxType == "Uniform") {
        const std::string part = ReadValue("particles:");
        filename += "_particle:" + part;
    }
    filename += "_" + GenSurface::DirectionTag();
    if (!Layout::Tag().empty()) filename += "_" + Layout::Tag();
    if (triggerType != "gamma") filename += "_trigger-" + triggerType;
    filename += ".txt";
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
        outDir += "_" + GenSurface::DirectionTag();
        if (!Layout::Tag().empty()) outDir += "_" + Layout::Tag();
        if (triggerType != "gamma") outDir += "_trigger-" + triggerType;
        outDir = sanitize(outDir);
        PostProcessing postProcessing(outDir, part);

        postProcessing.ExtractNtData();
        if (Emin < Emax) {
            if (genSurface.IsIsotropic())
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
