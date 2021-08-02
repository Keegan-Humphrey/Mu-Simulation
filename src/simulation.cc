
#include <G4MTRunManager.hh>
#include <FTFP_BERT.hh>
#include <G4StepLimiterPhysics.hh>
#include <G4UIExecutive.hh>
#include <G4VisExecutive.hh>
#include <tls.hh>

#include "action.hh"
#include "geometry/Construction.hh"
#include "geometry/Earth.hh"
#include "physics/Units.hh"
#include "ui.hh"
#include "PhysicsList.hh"
#include "MuonDataController.hh"

#include "G4GenericBiasingPhysics.hh"

#include "util/command_line_parser.hh"
#include "util/error.hh"

//__Main Function: Simulation___________________________________________________________________
int main(int argc, char* argv[]) {
  using namespace MATHUSLA;
  using namespace MATHUSLA::MU;

  using util::cli::option;

  option help_opt    ('h', "help",     "MATHUSLA Muon Simulation",  option::no_arguments);
  option gen_opt     ('g', "gen",      "Generator",                 option::required_arguments);
  option det_opt     ('d', "det",      "Detector",                  option::required_arguments);
  option shift_opt   (0,   "shift",    "Shift Last Earth Layer",    option::required_arguments);
  option data_opt    ('o' ,"out",      "Data Output Directory",     option::required_arguments);
  option export_opt  ('E', "export",   "Export Output Directory",   option::required_arguments);
  option script_opt  ('s', "script",   "Custom Script",             option::required_arguments);
  option events_opt  ('e', "events",   "Event Count",               option::required_arguments);
  option save_all_opt(0,   "save_all", "Save All Generator Events", option::no_arguments);
  option cut_save_opt(0,   "cut_save", "Save Events With Digi Cuts",option::no_arguments);
  option bias_opt    (0,   "bias",     "Bias Muon Nuclear Interactions in Earth Volume", option::no_arguments);
  option five_body_muon_decay_opt('f', "five_muon", "Make 3-body muon decay 5-body",     option::no_arguments);
  option non_random_muon_decay_opt('n',"non_random", "Make 5-body muon decays in order", option::no_arguments);
  option vis_opt     ('v', "vis",      "Visualization",             option::no_arguments);
  option quiet_opt   ('q', "quiet",    "Quiet Mode",                option::no_arguments);
  option thread_opt  ('j', "threads",  "Multi-Threading Mode: Specify Optional number of threads (default: 2)", option::optional_arguments);

  //TODO: pass quiet argument to builder and action initiaization to improve quietness

  const auto script_argc = -1 + util::cli::parse(argv,
    {&help_opt, &gen_opt, &det_opt, &shift_opt, &data_opt, &export_opt, &script_opt, &events_opt,
     &save_all_opt, &cut_save_opt, &bias_opt, &five_body_muon_decay_opt, &non_random_muon_decay_opt, &vis_opt, &quiet_opt, &thread_opt});


  util::error::exit_when(script_argc && !script_opt.argument,
    "[FATAL ERROR] Illegal Forwarding Arguments:\n"
    "              Passing arguments to simulation without script is disallowed.\n");

  G4UIExecutive* ui = nullptr;
  if (argc == 1 || vis_opt.count) {
    ui = new G4UIExecutive(argc, argv);
    vis_opt.count = 1;
  }

  util::error::exit_when(script_opt.argument && events_opt.argument,
    "[FATAL ERROR] Incompatible Arguments:\n",
    "              A script OR an event count can be provided, but not both.\n");

  G4Random::setTheEngine(new CLHEP::RanecuEngine);
  G4Random::setTheSeed(time(nullptr));


  if (thread_opt.argument) {
    auto opt = std::string(thread_opt.argument);
    if (opt == "on") {
      thread_opt.count = 2;
    } else if (opt == "off" || opt == "0") {
      thread_opt.count = 1;
    } else {
      try {
        thread_opt.count = std::stoi(opt);
      } catch (...) {
        thread_opt.count = 2;
      }
    }
  } else if (!thread_opt.count) {
    thread_opt.count = 2;
  }
  auto run = new G4MTRunManager;
  thread_opt.count=1;
  std::cout << "Warning!!!!! You can only run one thread.  This doesn't work, otherwise." << std::endl;
  run->SetNumberOfThreads(thread_opt.count);
  std::cout << "Running " << thread_opt.count
            << (thread_opt.count > 1 ? " Threads" : " Thread") << "\n";

  run->SetPrintProgress(1000);
  run->SetRandomNumberStore(false);

  Units::Define();

  if (shift_opt.argument)
    Earth::LastShift(std::stold(shift_opt.argument) * m);

  G4bool fiveBodyMuonDecays = five_body_muon_decay_opt.count;
  G4bool randomize = !(non_random_muon_decay_opt.count);

  MuonDataController* controller = new MuonDataController();
  controller->setRandom(randomize);
  controller->setOn(fiveBodyMuonDecays);

  G4GenericBiasingPhysics* biasingPhysics = new G4GenericBiasingPhysics();
  biasingPhysics->Bias( "mu+" );
  biasingPhysics->Bias( "mu-" );

  if(fiveBodyMuonDecays){
    auto physics = new PhysicsList();
    run->SetUserInitialization(physics);
  } else if (bias_opt.count){
    auto physics = new FTFP_BERT;
    physics->RegisterPhysics( biasingPhysics );
    physics->RegisterPhysics(new G4StepLimiterPhysics);
    run->SetUserInitialization(physics);
  } else{
    auto physics = new FTFP_BERT;
    physics->RegisterPhysics(new G4StepLimiterPhysics);
    run->SetUserInitialization(physics);
    util::error::exit_when(!randomize,"You have set the flag -n so that the order of five-body muon decays are not random, but you have not set -f to turn on five-body muon decays. \n Turn on five-body muon decays and try again, or do not use the flag -n");
  }

  const auto detector = det_opt.argument ? det_opt.argument : "Box";
  const auto export_dir = export_opt.argument ? export_opt.argument : "";
  run->SetUserInitialization(new Construction::Builder(detector, export_dir, save_all_opt.count, cut_save_opt.count));

  const auto generator = gen_opt.argument ? gen_opt.argument : "basic";
  const auto data_dir = data_opt.argument ? data_opt.argument : "data";
  run->SetUserInitialization(new ActionInitialization(generator, data_dir));


  auto vis = new G4VisExecutive("Quiet");
  vis->Initialize();

  Command::Execute("/run/initialize",
                   "/control/saveHistory scripts/G4History",
                   "/control/stopSavingHistory");

  Command::Execute(quiet_opt.count ? "/control/execute scripts/settings/quiet"
                                   : "/control/execute scripts/settings/verbose");

  if (vis_opt.count) {
    Command::Execute("/control/execute scripts/settings/init_vis");
    if (ui->IsGUI())
      Command::Execute("/control/execute scripts/settings/init_gui");
  }

  if (script_opt.argument) {
    util::error::exit_when(script_argc % 2,
      "[FATAL ERROR] Illegal Number of Script Forwarding Arguments:\n",
      "              Inputed ", script_argc, " arguments but forward arguments must be key-value pairs.\n");

    const auto script_path = std::string(script_opt.argument);
    if (script_argc) {
      for (std::size_t i{}; i < script_argc; i += 2) {
        Command::Execute("/control/alias " + std::string(argv[i + 1]) + " " + std::string(argv[i + 2]));
      }
      Command::Execute("/control/execute " + script_path);
    } else {
      Command::Execute("/control/execute " + script_path);
    }
  } else if (events_opt.argument) {
    Command::Execute("/run/beamOn " + std::string(events_opt.argument));
  }

  if (ui) {
    ui->SessionStart();
    delete ui;
  }

  delete vis;
  delete run;
  return 0;
}
//----------------------------------------------------------------------------------------------
