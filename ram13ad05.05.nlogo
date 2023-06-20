;; Simulation of AMR in Senegal poultry production
; Bacterial growth model based on Graesboll et al 2014
; Immunity model based on Mochan et al 2014

; CONFIGURE VARIABLES

; define global variables - strain characteristics
; 6 strains Ec EColi Sa Salmonella Au Autre/Other, each in both NR Non-resistant/susceptible and R Resistant forms
; each with standard variables (set in Interface):
; xm - max strain growth rate,
; ec50prim - concentration at which growth rate half of xm for primary antimicrobial (Tetracycline),
; gammaprim - Hill-coefficient for primary antimicrobial (Tetracycline),
; ec50seco - concentration at which growth rate half of xm for secondary antimicrobial (fluoroquinolones),
; gammaseco - Hill-coefficient for secondary antimicrobial (fluoroquinolones)
; start - starting quantities of strains (chick seeder populations)
; cap - carrying capacity of chicken's intestines
; fitness-cost - parameter establishing lower growth rate for R Resistant strains in absence of antimicrobials
; simulation-duration - determines length of production cycle (initial value 1008 hours = 42 days)
; stage - tracks stages in chicken production cycle
; stagephagorateadjust - scale of stage effect on phagorate (<1)
; Deadcounter - count of dead chickens collected
; Deadmass - total mass of chickens that die
; prim-tracker, seco-tracker - trackers for quantities of antimicrobials used over time
; patchtotalsPathoR-tracker - tracker for total patch pathogenic resistant bacteria over time
; response-counter - counter for responsive treatment strategies
; withdrawmoment - withdrawal period for antibiotics, set as hours before slaughter
; clean-counter - counts number of times CleanCage procedure executed
; clean-disposalPathoRsum - sums amount of pathogenic resistant bacteria disposed of over executions of CleanCage procedure
; carcass-disposalPathoRsum - sums amount of pathogenic resistant bacteria disposed of as dead birds removed
; season (configured in interface) - term denoting temperature / humidity affecting bacterial growth on patches
; feeders, not-feeders - agentsets of patches comprising feeders / not feeders
; initialpen - agentset of patches comprising initial rearing pen
; apparentEc, apparentSa, deathcauseEc, deathcauseSa - variables for monitoring apparent infections & cause of death
; AMUstance (configured in Interface) - variable for setting AMU strategy
; immunegrowthintersect - term dictating relationship of immune response to growth
; moveprob - p value for chicken moving in step, after starter phase
; mortcalib - calibration variable, observed mean final count chickens

Globals[EcNRxm EcNRec50prim EcNRgammaprim EcNRec50seco EcNRgammaseco
  EcRxm EcRec50prim EcRgammaprim EcRec50seco EcRgammaseco
  SaNRxm SaNRec50prim SaNRgammaprim SaNRec50seco SaNRgammaseco
  SaRxm SaRec50prim SaRgammaprim SaRec50seco SaRgammaseco
  AuNRxm AuNRec50prim AuNRgammaprim AuNRec50seco AuNRgammaseco
  AuRxm AuRec50prim AuRgammaprim AuRec50seco AuRgammaseco
  startEcNR startEcR startSaNR startSaR startAuNR startAuR
  cap
  simulation-duration stage stagephagorateadjust Deadcounter Deadmass
  prim-tracker seco-tracker patchtotalsPathoR-tracker response-counter withdrawmoment
  clean-counter clean-disposalPathoRsum carcass-disposalPathoRsum
  feeders not-feeders initialpen
  apparentEc apparentSa deathcauseEc deathcauseSa
  immunegrowthintersect
  mortcalib
]

; define chickens and variables
; quantEcNR, quantEcR, quantSaNR, quantSaR, quantAuNR, quantAuR - per chicken intestinal strain populations
; cprim, cseco - antimicrobial chicken concentrations for primary antimicrobial (Tetracycline) and secondary antimicrobial (fluoroquinolones)
; cutoffval - cutoff strain population value below which removal may occur, removprob - probability of removal of a strain when below cutoffval
; excreterate - excretion rate from intestine, intakerate - intake rate to intestine
; contamweight, surroundweight - weights accorded to infection probabilities
; contamscale, surroundscale - terms determining scale of infections
; surrounddelay - time for exposure to infection, determines length of surroundTracker list
; surroundTrackerEcNR, surroundTrackerEcR, surroundTrackerSaNR, surroundTrackerSaR - list terms for tracking exposure to surrounding infections
; infectprobs - overall term controlling how likely infections are to occur
; infectscales - overall term controlling scale of infections
; lungquantEcNR, lungquantEcR, lungquantSaNR, lungquantSaR - per chicken lungs pathogenic strain populations
; lunggrowthprop - proportion determining bacterial growth in lung relative to intestine
; lungmechclearrate, lungphagoclearrate -  rate at which the innate immune system initially controls the bacteria in the lung respectively through: mechanical clearance via cilia and mucus; phagocytosis by alveolar macrophages
; lungimmunesteady - steady-state value of lungquantImmune
; lungquantImmune - per chicken lung population of phagocyte immune cells
; Ecphagorate, Saphagorate - rate at which immune cells phagocytose pathogenic strain types, E Coli and Salmonella respectively
; phagoinhib - term for inhibition of phagocyte phagocytosis
; immuneinfluxrate - rate at which phagocyte immune cells initially influx to lungs
; immunedelay - time for immune cell response (time to activation of phagocytes), determines length of immuneTracker list
; immuneTracker - list term for tracking immune cell growth signal per tick
; epiDamage - damage to the epithelium of the lungs
; damageincreaserate, damagerepairrate - rates at which epiDamage: increase proportional to pathogenic strain populations; decrease through body repair
; home-feeder - feeder around which chicken moves
; mass - chicken weight
; maxmass - slaughter weight towards which chickens grow
; vax - vaccination status
; viva - death status
; countprim, countseco - treatment timers for primary and secondary antimicrobials

breed[chickens chicken]
chickens-own[quantEcNR quantEcR quantSaNR quantSaR quantAuNR quantAuR cprim cseco cutoffval removprob excreterate intakerate
  contamweight surroundweight contamscale surroundscale
  surrounddelay surroundTrackerEcNR surroundTrackerEcR surroundTrackerSaNR surroundTrackerSaR
  infectprobs infectscales
  lungquantEcNR lungquantEcR lungquantSaNR lungquantSaR lunggrowthprop lungmechclearrate lungphagoclearrate
  lungimmunesteady lungquantImmune Ecphagorate Saphagorate phagoinhib immuneinfluxrate immunedelay immuneTracker
  epiDamage damageincreaserate damagerepairrate
  home-feeder mass maxmass vax viva countprim countseco]

; define patch variables
; patchEcNR, patchEcR, patchSaNR, patchSaR, patchAuNR, patchAuR - per patch strain populations
; chickenCount - number of chickens on a patch
; patchgrowthprop - proportion determining bacterial growth on patch relative to intestines

patches-own[patchEcNR patchEcR patchSaNR patchSaR patchAuNR patchAuR chickenCount patchgrowthprop]

; CONFIGURE SETUP PROCEDURE
; -obs suffix used for variables set directly in Interface, otherwise -mean, -sd, -min, -max etc

to setup
  clear-all

  ; set cage size - 1 patch sides 16th m2 0.16667*0.16667cm (0.02778 m2 - 1/36th m2)
  ; standard cage 10m width, length determined by number birds at 10 per m2 (=> 20m length for 2000 birds)
  ; simulation for 500 birds / quarter standard size

  resize-world -30 29 -15 14

  __change-topology False False

  ; set global variables for all strains
  ; nb xm in R strains subject to fitness cost

  set EcNRxm random-normal EcNRxm_mean EcNRxm_sd
  set EcRxm (random-normal EcRxm_mean EcRxm_sd) * ((100 - random (EcRxm_maxfitcost% + 1)) / 100)
  set SaNRxm random-normal SaNRxm_mean SaNRxm_sd
  set SaRxm (random-normal SaRxm_mean SaRxm_sd) * ((100 - random (SaRxm_maxfitcost% + 1)) / 100)
  set AuNRxm random-normal AuNRxm_mean AuNRxm_sd
  set AuRxm (random-normal AuRxm_mean AuRxm_sd) * ((100 - random (AuRxm_maxfitcost% + 1)) / 100)

  set EcNRec50prim (random-float (EcNRec50primmax - EcNRec50primmin)) + EcNRec50primmin
  set EcRec50prim (random-float (EcRec50primmax - EcRec50primmin)) + EcRec50primmin
  set SaNRec50prim (random-float (SaNRec50primmax - SaNRec50primmin)) + SaNRec50primmin
  set SaRec50prim (random-float (SaRec50primmax - SaRec50primmin)) + SaRec50primmin
  set AuNRec50prim (random-float (AuNRec50primmax - AuNRec50primmin)) + AuNRec50primmin
  set AuRec50prim (random-float (AuRec50primmax - AuRec50primmin)) + AuRec50primmin

  set EcNRgammaprim (random-float (EcNRgammaprimmax - EcNRgammaprimmin)) + EcNRgammaprimmin
  set EcRgammaprim (random-float (EcRgammaprimmax - EcRgammaprimmin)) + EcRgammaprimmin
  set SaNRgammaprim (random-float (SaNRgammaprimmax - SaNRgammaprimmin)) + SaNRgammaprimmin
  set SaRgammaprim (random-float (SaRgammaprimmax - SaRgammaprimmin)) + SaRgammaprimmin
  set AuNRgammaprim (random-float (AuNRgammaprimmax - AuNRgammaprimmin)) + AuNRgammaprimmin
  set AuRgammaprim (random-float (AuRgammaprimmax - AuRgammaprimmin)) + AuRgammaprimmin

  set EcNRec50seco (random-float (EcNRec50secomax - EcNRec50secomin)) + EcNRec50secomin
  set EcRec50seco (random-float (EcRec50secomax - EcRec50secomin)) + EcRec50secomin
  set SaNRec50seco (random-float (SaNRec50secomax - SaNRec50secomin)) + SaNRec50secomin
  set SaRec50seco (random-float (SaRec50secomax - SaRec50secomin)) + SaRec50secomin
  set AuNRec50seco (random-float (AuNRec50secomax - AuNRec50secomin)) + AuNRec50secomin
  set AuRec50seco (random-float (AuRec50secomax - AuRec50secomin)) + AuRec50secomin

  set EcNRgammaseco (random-float (EcNRgammasecomax - EcNRgammasecomin)) + EcNRgammasecomin
  set EcRgammaseco (random-float (EcRgammasecomax - EcRgammasecomin)) + EcRgammasecomin
  set SaNRgammaseco (random-float (SaNRgammasecomax - SaNRgammasecomin)) + SaNRgammasecomin
  set SaRgammaseco (random-float (SaRgammasecomax - SaRgammasecomin)) + SaRgammasecomin
  set AuNRgammaseco (random-float (AuNRgammasecomax - AuNRgammasecomin)) + AuNRgammasecomin
  set AuRgammaseco (random-float (AuRgammasecomax - AuRgammasecomin)) + AuRgammasecomin

  ; set intestinal capacity of chickens

  set cap 1.0E9

  ; set duration, stage, trackers etc.

  ; duration 6 weeks, 42 days - 1,008 hours

  set simulation-duration 1008

  set stage "starter"
  set stagephagorateadjust stagephagorateadjustobs
  set Deadcounter 0
  set Deadmass 0

  set apparentEc 0
  set apparentSa 0
  set deathcauseEc 0
  set deathcauseSa 0

  ; trackers duration + 1 week for delayed slaughter

  set prim-tracker (n-values (simulation-duration + 168) [0])
  set seco-tracker (n-values (simulation-duration + 168) [0])
  set patchtotalsPathoR-tracker (n-values (simulation-duration + 168) [0])

  set response-counter 0
  set withdrawmoment (7 * 24) ; 7 days
  set clean-counter 0
  set clean-disposalPathoRsum 0
  set carcass-disposalPathoRsum 0

  ; define feeders location (1 per 50 birds, approximately evenly distributed)

  set feeders (patch-set patch -20 -5 patch -20 5
    patch -10 -5 patch -10 5
    patch 0 -5 patch 0 5
    patch 10 -5 patch 10 5
    patch 20 -5 patch 20 5)

  set not-feeders patches with [not member? self feeders]

  set immunegrowthintersect immunegrowthintersectobs

  ; define starter circles (https://elevage.sec.gouv.sn/sites/default/files/GUIDE%20DELABORATION%20DE%20PROJET%20POULETS%20DE%20CHAIR%20VF.pdf)

  InitialPenSet

  ; set calibration variables

  set mortcalib 50

  ; deploy chickens - approx 10 chicken per m2 determines number of chickens / chickens-to-cage size ratio
  ; chickens initially kept within smaller rearing area - see MoveStep procedure

  create-chickens 500
  ask chickens [move-to one-of initialpen]

  ; set chickens general initial values
  ; nb cutoffval set to midrange 50,000, excreterate set to midrange 0.01, intakerate set to midrange 0.001
  ; infectprobs - value arbitrary, designed to help calibrate model to observed disease incidence
  ; infectscales - set with reference to initial inoculation values used in Mochan et al 2014 (1.0E5 and 1.0E6)
  ; nb lunggrowthprop initial 0.5
  ; nb lungmechclearrate, lungphagoclearrate initial values taken in range proposed by Mochan et al 2014 (respectively 1.0E4 - 1.0E9, 1.0E5 - 1.0E7)
  ; nb lungimmunesteady initial value in midrange proposed by Mochan et al 2014 (1.0E2 - 1.0E4)
  ; nb Ecphagorate, Saphagorate initial values taken from within range proposed by Mochan et al 2014 (1.0E-9 - 1.0E-4)
  ; nb phagoinhib initial value in range midrange proposed by Mochan et al 2014 (1.0E-13 - 1.0E-11)
  ; nb immuneinfluxrate initial value in range proposed by Mochan et al 2014 (1.0E5 - 1.0E8)
  ; nb immunedelay initial value in range proposed by Mochan et al 2014 (0.01 - 10)
  ; nb damageincreaserate initial values taken from within range proposed by Mochan et al 2014 (1.0E-9 - 1.0E-7)
  ; nb damagerepairrate initial values taken from within range proposed by Mochan et al 2014 (0.01–1)
  ; nb starting mass approx. 200g from Mouffok et al 2019
  ; nb maxmass - 2000 - 2kg

  ask chickens [set shape "triangle" set color yellow set size 1 set cprim 0 set cseco 0 set cutoffval cutoffvalobs set removprob removprobobs set excreterate excreterateobs set intakerate intakerateobs
    set contamweight contamweightobs set surroundweight surroundweightobs set contamscale contamscaleobs set surroundscale surroundscaleobs
    set surrounddelay surrounddelayobs set surroundTrackerEcNR (n-values surrounddelay [0]) set surroundTrackerEcR (n-values surrounddelay [0]) set surroundTrackerSaNR (n-values surrounddelay [0]) set surroundTrackerSaR (n-values surrounddelay [0])
    set infectprobs infectprobsobs set infectscales infectscalesobs
    set lungquantEcNR 0 set lungquantEcR 0 set lungquantSaNR 0 set lungquantSaR 0 set lunggrowthprop lunggrowthpropobs set lungmechclearrate lungmechclearrateobs set lungphagoclearrate lungphagoclearrateobs
    set lungimmunesteady lungimmunesteadyobs set lungquantImmune (lungimmunesteady) set Ecphagorate Ecphagorateobs set Saphagorate Saphagorateobs
    set phagoinhib phagoinhibobs set immuneinfluxrate immuneinfluxrateobs set immunedelay immunedelayobs set immuneTracker (n-values immunedelay [0])
    set epiDamage 0 set damageincreaserate damageincreaserateobs set damagerepairrate damagerepairrateobs
    set home-feeder one-of feeders set mass 200 set maxmass 2000 set vax "No" set viva "Alive"]

  ; set chickens individual values
  ; presence/varying quantities of different strains in each chicken
  ; nb only varying proportions of NR / R Ec and Sa strains for now
  ; nb proposing starting point at 0.25 of capacity occupied

  ask chickens [
    set quantEcR propEc * startcap * ((random-float (maxRinEc - minRinEc)) + minRinEc)
    set quantEcNR (propEc * startcap) - quantEcR
    set quantSaR propSa * startcap * ((random-float (maxRinSa - minRinSa)) + minRinSa)
    set quantSaNR (propSa * startcap) - quantSaR
    set quantAuR (1 - (propEc + propSa)) * startcap * RinAu
    set quantAuNR (1 - (propEc + propSa)) * startcap * (1 - RinAu)
  ]

  ; set patches initial values
  ; patchgrowthprop - arbitrary
  ; season - set in interface

  ask patches [set pcolor 65 set patchEcNR 0 set patchEcR 0 set patchSaNR 0 set patchSaR 0 set patchAuNR 0 set patchAuR 0 set chickenCount (count turtles-here) set patchgrowthprop patchgrowthpropobs]

  ; show feeder patches

  ask feeders [set pcolor brown]

  reset-ticks

end

; CONFIFURE GO PROCEDURE

to go

  ; end procedure
  ; nb tick is 1 hour
  ; nb mortcalib for calibration metrics

  ; stop if zero chickens

  if (count chickens = 0) [

    set stage "all dead"

    ; handle any unremoved dead birds

    ask chickens [if (viva = "Dead") [DisposeDead]]

    ; report outcomes

    set clean-disposalPathoRsum (clean-disposalPathoRsum + sum [patchEcR] of patches + sum [patchSaR] of patches)

    set mortcalib ((count chickens) - 450) * ((count chickens) - 450)

    carefully [Report-Outcomes] [output-print error-message]

    stop
  ]

  ; stop at term provided that 90% of birds at target weight

  if ticks >= simulation-duration and ((count chickens with [mass > maxmass]) / (count chickens) > 0.9) [

    set stage "slaughter"

    ; handle any unremoved dead birds

    ask chickens [if (viva = "Dead") [DisposeDead]]

    ; report outcomes

    set clean-disposalPathoRsum (clean-disposalPathoRsum + sum [patchEcR] of patches + sum [patchSaR] of patches)

    set mortcalib ((count chickens) - 450) * ((count chickens) - 450)

    carefully [Report-Outcomes] [output-print error-message]

    stop
  ]

  ; hard stop if 1 week beyond term

  if ticks >= (simulation-duration + 168) [

    set stage "slaughter"

    ; handle any unremoved dead birds

    ask chickens [if (viva = "Dead") [DisposeDead]]

    ; report outcomes

    set clean-disposalPathoRsum (clean-disposalPathoRsum + sum [patchEcR] of patches + sum [patchSaR] of patches)

    set mortcalib ((count chickens) - 450) * ((count chickens) - 450)

    carefully [Report-Outcomes] [output-print error-message]

    stop
  ]

  ; adjust stage global variable (chicken life stage) - starter up to 14 days, then grower to 28 days then finisher to slaughter (see end above)

  if ticks >= 336 and ticks < 672 [set stage "grower"]

  if ticks >= 672 [set stage "finisher"]

  ; bacterial growth calculations - see subordinate procedure

  ask chickens [carefully [StrainsGrow] [print error-message]]

  ; intake - provisional, see subordinate procedure

  ask chickens [StrainsIntake]

  ; excrete - provisional, see subordinate procedure

  ask chickens [StrainsExcrete]

  ; remove strains

  ask chickens [if (quantEcNR < cutoffval) and ((random 11 / 10) < removprob) [set quantEcNR 0]
    if (quantEcR < cutoffval) and ((random 11 / 10) < removprob) [set quantEcR 0]
    if (quantSaNR < cutoffval) and ((random 11 / 10) < removprob) [set quantSaNR 0]
    if (quantSaR < cutoffval) and ((random 11 / 10) < removprob) [set quantSaR 0]
    if (quantAuNR < cutoffval) and ((random 11 / 10) < removprob) [set quantAuNR 0]
    if (quantAuR < cutoffval) and ((random 11 / 10) < removprob) [set quantAuR 0]]

  ; infection-related procedures

  ask chickens [if (viva = "Alive") [InfectionMech]]

  ; immunity-related procedures

  ask chickens [if (viva = "Alive") [ImmunityStep]]

  ; patch strain growth

  ask patches [carefully [PatchStrainGrow] [print error-message]]

  ; update patch colours - provisional, based on hygiene (total Ec & Sa excreted vs 0.5 * cap) green to black if any patch strain population * 0.0015 is less than cutoffval, light to dark red if above or equal

  ask not-feeders [UpdateColors]

  ; illness signs

  ask chickens [

    if (epiDamage > 0.2) and (viva = "Alive") [set color 15]

    if (epiDamage <= 0.2) and (viva = "Alive") [set color yellow]

  ]

  ; mortality procedures - death if too much lung damage, dead birds removed every 12 ticks

  ask chickens [

    if (viva = "Alive") [if (epiDamage > 0.8) [set color 0 set viva "Dead" set Deadmass (Deadmass + mass) set mass 0]]

    if (viva = "Alive") [if (mass < 200) [set color 0 set viva "Dead" set Deadmass (Deadmass + mass) set mass 0]]

  ]

  ask chickens [

    if (viva = "Dead") and ((ticks mod 12) = 0) [DisposeDead]

  ]

  ; grow weight - PROVISIONAL/NEEDS TESTING, approx. sigmoid with impact disease after incubation - see subordinate function

  ask chickens [if (viva = "Alive") [GrowMass]]

  ; move (including feeding every 12 hours) and update patch counts

  MoveStep

  ; treatment - provisional, see subordinate procedures (called using Interface)

  ; set treatment counters plus antimicrobial concentrations affecting relative bacterial growth rates

  ask chickens [if countprim > 0 [set countprim (countprim - 1)] if countseco > 0 [set countseco (countseco - 1)]]

  if AMUstance != "Responsive quick stop" and AMUstance != "Prophylactic" [

    ask chickens [if countprim = 0 [set cprim 0] if countseco = 0 [set cseco 0]]

  ]

  ; cleaning enclosure - provisional, see subordinate procedures (called using Interface)

  ; INTERVENTIONS (called from AMUstance in Interface)

  if AMUstance = "Responsive quick stop" [

    if (count (chickens with [epiDamage > 0.2]) >= 3) and response-counter < 48 [

      ask chickens [set cprim 10]

      set response-counter (response-counter + 1)

    ]

    if (count (chickens with [epiDamage > 0.2]) >= 3) and response-counter >= 48 and response-counter < 96 [

      ask chickens [set cprim 40]

      set response-counter (response-counter + 1)

    ]

    if (count (chickens with [epiDamage > 0.2]) >= 3) and response-counter >= 96 [

      ask chickens [set cprim 40 set cseco 40]

      set response-counter (response-counter + 1)

    ]

    if (count (chickens with [epiDamage > 0.2]) < 3) [

      ask chickens [set cprim 0 set cseco 0]

      set response-counter 0

    ]

  ]

  if AMUstance = "Responsive complete" [

    if (count (chickens with [epiDamage > 0.2]) >= 3) and (response-counter = 0 or (response-counter = "second" and (sum [cseco] of chickens) = 0)) [

      ask chickens [set cprim 10 set countprim 168]

      set response-counter "first"

    ]

    if (count (chickens with [epiDamage > 0.2]) >= 3) and response-counter = "first" and (sum [cprim] of chickens) = 0 [

      ask chickens [set cseco 40 set countseco 72]

      set response-counter "second"

    ]

    if ticks > (simulation-duration - withdrawmoment) [

      ask chickens [set cprim 0 set cseco 0]

      set response-counter "withdrawn"

    ]

  ]

  if AMUstance = "Prophylactic" [

    if response-counter != "withdrawn" [ask chickens [set cprim 1]]

    if ticks > (simulation-duration - withdrawmoment) [ask chickens [set cprim 0] set response-counter "withdrawn"]

  ]

  if AMUstance = "Clean responsive" [

    if (count (chickens with [epiDamage > 0.2]) >= 3) and ticks mod 72 = 0 [

      CleanCageAuto

    ]

  ]

  if AMUstance = "Regular clean" [

    if ticks mod 168 = 0 [

      CleanCageAuto

    ]

  ]

  if AMUstance = "Responsive quick stop + clean" [

    if (count (chickens with [epiDamage > 0.2]) >= 3) and response-counter < 48 [

      ask chickens [set cprim 10]

      set response-counter (response-counter + 1)

    ]

    if (count (chickens with [epiDamage > 0.2]) >= 3) and response-counter >= 48 and response-counter < 96 [

      ask chickens [set cprim 40]

      set response-counter (response-counter + 1)

    ]

    if (count (chickens with [epiDamage > 0.2]) >= 3) and response-counter >= 96 [

      ask chickens [set cprim 40 set cseco 40]

      set response-counter (response-counter + 1)

    ]

    if (count (chickens with [epiDamage > 0.2]) < 3) [

      ask chickens [set cprim 0 set cseco 0]

      set response-counter 0

    ]

    if (count (chickens with [epiDamage > 0.2]) >= 3) and ticks mod 72 = 0 [

      CleanCageAuto

    ]

  ]

  if AMUstance = "Responsive complete + reg. clean" [

    if (count (chickens with [epiDamage > 0.2]) >= 3) and (response-counter = 0 or (response-counter = "second" and (sum [cseco] of chickens) = 0)) [

      ask chickens [set cprim 10 set countprim 168]

      set response-counter "first"

    ]

    if (count (chickens with [epiDamage > 0.2]) >= 3) and response-counter = "first" and (sum [cprim] of chickens) = 0 [

      ask chickens [set cseco 40 set countseco 72]

      set response-counter "second"

    ]

    if ticks > (simulation-duration - withdrawmoment) [

      ask chickens [set cprim 0 set cseco 0]

      set response-counter "withdrawn"

    ]

    if ticks mod 168 = 0 [

      CleanCageAuto

    ]

  ]

  ; update trackers & monitors

  set prim-tracker (replace-item ticks prim-tracker sum [cprim] of chickens)

  set seco-tracker (replace-item ticks seco-tracker sum [cseco] of chickens)

  set patchtotalsPathoR-tracker (replace-item ticks patchtotalsPathoR-tracker (sum [patchEcR] of patches + sum [patchSaR] of patches))

  set apparentEc count chickens with [epiDamage > 0.2 and ((lungquantEcNR + lungquantEcR) > (lungquantSaNR + lungquantSaR))]
  set apparentSa count chickens with [epiDamage > 0.2 and ((lungquantSaNR + lungquantSaR) > (lungquantEcNR + lungquantEcR))]

  tick

end



; CONFIGURE SUBORDINATE PROCEDURES

; define initial pen areas

to InitialPenSet

  let demarrecercle1 [patches in-radius 12] of patch 0 0

  set initialpen (patch-set demarrecercle1)

end

; define procedure for bacterial growth

to StrainsGrow

  ; per strain growth rate Hij given antimicrobial concentrations
  ; nb sequential terms deal with effect of primary (Tetracycline) then secondary (fluoroquinolones)

  let HijEcNR (EcNRxm * (1 - ((cprim ^ EcNRgammaprim) / ((EcNRec50prim ^ EcNRgammaprim) + (cprim ^ EcNRgammaprim)))) * (1 - ((cseco ^ EcNRgammaseco) / ((EcNRec50seco ^ EcNRgammaseco) + (cseco ^ EcNRgammaseco)))))
  let HijEcR (EcRxm * (1 - ((cprim ^ EcRgammaprim) / ((EcRec50prim ^ EcRgammaprim) + (cprim ^ EcRgammaprim)))) * (1 - ((cseco ^ EcRgammaseco) / ((EcRec50seco ^ EcRgammaseco) + (cseco ^ EcRgammaseco)))))
  let HijSaNR (SaNRxm * (1 - ((cprim ^ SaNRgammaprim) / ((SaNRec50prim ^ SaNRgammaprim) + (cprim ^ SaNRgammaprim)))) * (1 - ((cseco ^ SaNRgammaseco) / ((SaNRec50seco ^ SaNRgammaseco) + (cseco ^ SaNRgammaseco)))))
  let HijSaR (SaRxm * (1 - ((cprim ^ SaRgammaprim) / ((SaRec50prim ^ SaRgammaprim) + (cprim ^ SaRgammaprim)))) * (1 - ((cseco ^ SaRgammaseco) / ((SaRec50seco ^ SaRgammaseco) + (cseco ^ SaRgammaseco)))))
  let HijAuNR (AuNRxm * (1 - ((cprim ^ AuNRgammaprim) / ((AuNRec50prim ^ AuNRgammaprim) + (cprim ^ AuNRgammaprim)))) * (1 - ((cseco ^ AuNRgammaseco) / ((AuNRec50seco ^ AuNRgammaseco) + (cseco ^ AuNRgammaseco)))))
  let HijAuR (AuRxm * (1 - ((cprim ^ AuRgammaprim) / ((AuRec50prim ^ AuRgammaprim) + (cprim ^ AuRgammaprim)))) * (1 - ((cseco ^ AuRgammaseco) / ((AuRec50seco ^ AuRgammaseco) + (cseco ^ AuRgammaseco)))))

  ; sum existing strain counts

  let sumStrains (quantEcNR + quantEcR + quantSaNR + quantSaR + quantAuNR + quantAuR)

  ; per strain growth Gij given competitive positions
  ; nb query appropriate implementation of second order differential equations? (doesn't the sum term tend to zero here?)

  let GijEcNR (HijEcNR * quantEcNR * (((cap - quantEcNR) * (cap - sumStrains)) / (cap ^ 2)))
  let GijEcR (HijEcR * quantEcR * (((cap - quantEcR) * (cap - sumStrains)) / (cap ^ 2)))
  let GijSaNR (HijSaNR * quantSaNR * (((cap - quantSaNR) * (cap - sumStrains)) / (cap ^ 2)))
  let GijSaR (HijSaR * quantSaR * (((cap - quantSaR) * (cap - sumStrains)) / (cap ^ 2)))
  let GijAuNR (HijAuNR * quantAuNR * (((cap - quantAuNR) * (cap - sumStrains)) / (cap ^ 2)))
  let GijAuR (HijAuR * quantAuR * (((cap - quantAuR) * (cap - sumStrains)) / (cap ^ 2)))

  ; implement growth in chickens

  set quantEcNR (quantEcNR + GijEcNR)
  set quantEcR (quantEcR + GijEcR)
  set quantSaNR (quantSaNR + GijSaNR)
  set quantSaR (quantSaR + GijSaR)
  set quantAuNR (quantAuNR + GijAuNR)
  set quantAuR (quantAuR + GijAuR)

end

; define procedure for intake from patch

to StrainsIntake

  ; calculate intake amounts per strain

  let intakeEcNR (patchEcNR * intakerate)
  let intakeEcR (patchEcR * intakerate)
  let intakeSaNR (patchSaNR * intakerate)
  let intakeSaR (patchSaR * intakerate)
  let intakeAuNR (patchAuNR * intakerate)
  let intakeAuR (patchAuR * intakerate)

  ; add to chicken strain counts (dividing among chickens on patch)

  set quantEcNR (quantEcNR + (intakeEcNR / chickenCount))
  set quantEcR (quantEcR + (intakeEcR / chickenCount))
  set quantSaNR (quantSaNR + (intakeSaNR / chickenCount))
  set quantSaR (quantSaR + (intakeSaR / chickenCount))
  set quantAuNR (quantAuNR + (intakeAuNR / chickenCount))
  set quantAuR (quantAuR + (intakeAuR / chickenCount))

  ; subtract from patch strain counts where chickens stand

  if chickenCount > 0 [
     set patchEcNR (patchEcNR - intakeEcNR)
     set patchEcR (patchEcR - intakeEcR)
     set patchSaNR (patchSaNR - intakeSaNR)
     set patchEcR (patchSaR - intakeSaR)
     set patchAuNR (patchAuNR - intakeAuNR)
     set patchAuR (patchAuR - intakeAuR)
  ]

end

; define procedure for intake excrete to patch

to StrainsExcrete

  ; calculate excrete amounts per strain

  let excreteEcNR (quantEcNR * excreterate)
  let excreteEcR (quantEcR * excreterate)
  let excreteSaNR (quantSaNR * excreterate)
  let excreteSaR (quantSaR * excreterate)
  let excreteAuNR (quantAuNR * excreterate)
  let excreteAuR (quantAuR * excreterate)

  ; subtract from chicken strain counts

  set quantEcNR (quantEcNR - excreteEcNR)
  set quantEcR (quantEcR  - excreteEcR)
  set quantSaNR (quantSaNR  - excreteSaNR)
  set quantSaR (quantSaR  - excreteSaR)
  set quantAuNR (quantAuNR - excreteAuNR)
  set quantAuR (quantAuR - excreteAuR)

  ; add to patch strain counts

  set patchEcNR (patchEcNR + excreteEcNR)
  set patchEcR (patchEcR + excreteEcR)
  set patchSaNR (patchSaNR + excreteSaNR)
  set patchSaR (patchSaR + excreteSaR)
  set patchAuNR (patchAuNR + excreteAuNR)
  set patchAuR (patchAuR + excreteAuR)

end

; define infection-related procedures

to InfectionMech

  ;; determine probability of lung infection as lab-style dose (loosely inspired 1. by Becker et al 2022 dual vectors - seeder bird / contaminated pen and 2. papers on bird-to-bird & aerosol transmission)
  ;; function here being to represent infections, relating bacterial growth and disease submodels to each other in a manner capable of being calibrated to observed morbidity, mortality rates
  ;; PROVISIONAL/NEEDS TESTING

  ; contamination probability and scale determined as a function of patch bacteria

  let patchMax (cap * excreterate * ((count chickens) + Deadcounter))

  let probEcNRcontam 0
  let probEcRcontam 0
  let probSaNRcontam 0
  let probSaRcontam 0

  if (patchEcNR > infectscales) [set probEcNRcontam (patchEcNR / patchMax)]
  if (patchEcR > infectscales) [set probEcRcontam (patchEcR / patchMax)]
  if (patchSaNR > infectscales) [set probSaNRcontam (patchSaNR / patchMax)]
  if (patchSaR > infectscales) [set probSaRcontam (patchSaR / patchMax)]

  let scaleEcNRcontam (patchEcNR / patchMax)
  let scaleEcRcontam (patchEcR / patchMax)
  let scaleSaNRcontam (patchSaNR / patchMax)
  let scaleSaRcontam (patchSaR / patchMax)

  ;if who = 0 [ask chicken 0 [

    ;print (word "contam probs" probEcNRcontam probEcRcontam probSaNRcontam probSaRcontam)

    ;print (word "contam scales" scaleEcNRcontam scaleEcRcontam scaleSaNRcontam scaleSaRcontam)

  ;]]

  ; surroundings probability and scale determined as a function of neighbouring or encountered birds already infected (modelling aerosol)

  let surroundingEcNR ((sum [lungquantEcNR] of turtles-on neighbors) + (sum [lungquantEcNR] of turtles-here) - lungquantEcNR)
  let surroundingEcR ((sum [lungquantEcR] of turtles-on neighbors) + (sum [lungquantEcR] of turtles-here) - lungquantEcR)
  let surroundingSaNR ((sum [lungquantSaNR] of turtles-on neighbors) + (sum [lungquantSaNR] of turtles-here) - lungquantSaNR)
  let surroundingSaR ((sum [lungquantSaR] of turtles-on neighbors) + (sum [lungquantSaR] of turtles-here) - lungquantSaR)

  set surroundTrackerEcNR fput surroundingEcNR but-last surroundTrackerEcNR
  set surroundTrackerEcR fput surroundingEcR but-last surroundTrackerEcR
  set surroundTrackerSaNR fput surroundingSaNR but-last surroundTrackerSaNR
  set surroundTrackerSaR fput surroundingSaR but-last surroundTrackerSaR

  let recentsurroundingEcNR (sum(surroundTrackerEcNR))
  let recentsurroundingEcR (sum(surroundTrackerEcR))
  let recentsurroundingSaNR (sum(surroundTrackerSaNR))
  let recentsurroundingSaR (sum(surroundTrackerSaR))

  let probEcNRsurround 0
  let probEcRsurround 0
  let probSaNRsurround 0
  let probSaRsurround 0

  if (recentsurroundingEcNR > infectscales) [set probEcNRsurround (recentsurroundingEcNR / patchMax)]
  if (recentsurroundingEcR > infectscales) [set probEcRsurround (recentsurroundingEcR / patchMax)]
  if (recentsurroundingSaNR > infectscales) [set probSaNRsurround (recentsurroundingSaNR / patchMax)]
  if (recentsurroundingSaR > infectscales) [set probSaRsurround (recentsurroundingSaR / patchMax)]

  let scaleEcNRsurround (recentsurroundingEcNR / patchMax)
  let scaleEcRsurround (recentsurroundingEcR / patchMax)
  let scaleSaNRsurround (recentsurroundingSaNR / patchMax)
  let scaleSaRsurround (recentsurroundingSaR / patchMax)

  ;if who = 0 and ((surroundingEcNR + surroundingEcR + surroundingSaNR + surroundingSaR) > 0) [ask chicken 0 [

    ;print (word "surround probs" probEcNRsurround probEcRsurround probSaNRsurround probSaRsurround)

    ;print (word "surrounding bacteria" surroundingEcNR surroundingEcR surroundingSaNR surroundingSaR)

    ;print (word "surround scales" scaleEcNRsurround scaleEcRsurround scaleSaNRsurround scaleSaRsurround)

  ;]]

  ; apply vector / overall probability weights

  let allprobweights (contamweight + surroundweight)

  let contamweightapplied ((contamweight / allprobweights) * infectprobs)
  let surroundweightapplied ((surroundweight / allprobweights) * infectprobs)

  set probEcNRcontam (probEcNRcontam * contamweightapplied)
  set probEcRcontam (probEcRcontam * contamweightapplied)
  set probSaNRcontam (probSaNRcontam * contamweightapplied)
  set probSaRcontam (probSaRcontam * contamweightapplied)

  set probEcNRsurround (probEcNRsurround * surroundweightapplied)
  set probEcRsurround (probEcRsurround * surroundweightapplied)
  set probSaNRsurround (probSaNRsurround * surroundweightapplied)
  set probSaRsurround (probSaRsurround * surroundweightapplied)

  ; apply vector / overall scales

  set scaleEcNRcontam (scaleEcNRcontam * infectscales * contamscale)
  set scaleEcRcontam (scaleEcRcontam * infectscales * contamscale)
  set scaleSaNRcontam (scaleSaNRcontam * infectscales * contamscale)
  set scaleSaRcontam (scaleSaRcontam * infectscales * contamscale)

  set scaleEcNRsurround (scaleEcNRsurround * infectscales * surroundscale)
  set scaleEcRsurround (scaleEcRsurround * infectscales * surroundscale)
  set scaleSaNRsurround (scaleSaNRsurround * infectscales * surroundscale)
  set scaleSaRsurround (scaleSaRsurround * infectscales * surroundscale)

  ; distribute scales normally within bounds - nb issue with spikes at boundary of range

  set scaleEcNRcontam median (list 0 (random-normal scaleEcNRcontam infectscales) patchEcNR)
  set scaleEcRcontam median (list 0 (random-normal scaleEcRcontam infectscales) patchEcR)
  set scaleSaNRcontam median (list 0 (random-normal scaleSaNRcontam infectscales) patchSaNR)
  set scaleSaRcontam median (list 0 (random-normal scaleSaRcontam infectscales) patchSaR)

  set scaleEcNRsurround median (list 0 (random-normal scaleEcNRsurround infectscales) infectscales)
  set scaleEcRsurround median (list 0 (random-normal scaleEcRsurround infectscales) infectscales)
  set scaleSaNRsurround median (list 0 (random-normal scaleSaNRsurround infectscales) infectscales)
  set scaleSaRsurround median (list 0 (random-normal scaleSaRsurround infectscales) infectscales)

  ; distribute scales normally within bounds - arbitrary, see reporter function below - discarded to avoid recursion & issues handling zero values https://stackoverflow.com/questions/20230685/netlogo-how-to-make-sure-a-variable-stays-in-a-defined-range

  ;if patchEcNR > 0 [set scaleEcNRcontam (random-normal-in-bounds scaleEcNRcontam infectscales 0 patchEcNR)]
  ;if patchEcR > 0 [set scaleEcRcontam (random-normal-in-bounds scaleEcRcontam infectscales 0 patchEcR)]
  ;if patchSaNR > 0 [set scaleSaNRcontam (random-normal-in-bounds scaleSaNRcontam infectscales 0 patchSaNR)]
  ;if patchSaR > 0 [set scaleSaRcontam (random-normal-in-bounds scaleSaRcontam infectscales 0 patchSaR)]

  ;set scaleEcNRsurround (random-normal-in-bounds scaleEcNRsurround infectscales 0 infectscales)
  ;set scaleEcRsurround (random-normal-in-bounds scaleEcRsurround infectscales 0 infectscales)
  ;set scaleSaNRsurround (random-normal-in-bounds scaleSaNRsurround infectscales 0 infectscales)
  ;set scaleSaRsurround (random-normal-in-bounds scaleSaRsurround infectscales 0 infectscales)


  ; implement infections - contamination, surroundings

  if random-float 1 < probEcNRcontam [set lungquantEcNR (lungquantEcNR + scaleEcNRcontam)
    ;print (word "Chicken " who " infected with EcNR via contam vector (probability " probEcNRcontam " scale " scaleEcNRcontam ")")
  ]
  if random-float 1 < probEcRcontam [set lungquantEcR (lungquantEcR + scaleEcRcontam)
    ;print (word "Chicken " who " infected with EcR via contam vector (probability " probEcRcontam " scale " scaleEcRcontam ")")
  ]
  if random-float 1 < probSaNRcontam [set lungquantSaNR (lungquantSaNR + scaleSaNRcontam)
    ;print (word "Chicken " who " infected with SaNR via contam vector (probability " probSaNRcontam " scale " scaleSaNRcontam ")")
  ]
  if random-float 1 < probSaRcontam [set lungquantSaR (lungquantSaR + scaleSaRcontam)
    ;print (word "Chicken " who " infected with SaR via contam vector (probability " probSaRcontam " scale " scaleSaRcontam ")")
  ]

  if random-float 1 < probEcNRsurround [set lungquantEcNR (lungquantEcNR + scaleEcNRsurround)
    ;print (word "Chicken " who " infected with EcNR via surround vector (probability " probEcNRsurround " scale " scaleEcNRsurround ")")
  ]
  if random-float 1 < probEcRsurround [set lungquantEcR (lungquantEcR + scaleEcRsurround)
    ;print (word "Chicken " who " infected with EcR via surround vector (probability " probEcRsurround " scale " scaleEcRsurround ")")
  ]
  if random-float 1 < probSaNRsurround [set lungquantSaNR (lungquantSaNR + scaleSaNRsurround)
    ;print (word "Chicken " who " infected with SaNR via surround vector (probability " probSaNRsurround " scale " scaleSaNRsurround ")")
  ]
  if random-float 1 < probSaRsurround [set lungquantSaR (lungquantSaR + scaleSaRsurround)
    ;print (word "Chicken " who " infected with SaR via surround vector (probability " probSaRsurround " scale " scaleSaRsurround ")")
  ]

  ; implement consequences of infections for patches' bacterial populations? - not currently implemented

end

; reporter function for normal distribution within bounds

;to-report random-normal-in-bounds [mid dev mmin mmax]
  ;let result random-normal mid dev
  ;if result < mmin or result > mmax
    ;[ report random-normal-in-bounds mid dev mmin mmax ]
  ;report result
;end

; define immunity-related procedures

to ImmunityStep

  ;; nb alternative to infection step is strengthen growth, weaken innate immune and let exposure build up to Mochan model limit over time, maybe with threshold for immune response?

  ;; grow lung pathogenic bacteria populations

  ; determine growth rates  using same function as for intestine growth (ie including modification for antimicrobials), but subject to lung growth rate penalty

  let lunggrowrateEcNR (lunggrowthprop * (EcNRxm * (1 - ((cprim ^ EcNRgammaprim) / ((EcNRec50prim ^ EcNRgammaprim) + (cprim ^ EcNRgammaprim)))) * (1 - ((cseco ^ EcNRgammaseco) / ((EcNRec50seco ^ EcNRgammaseco) + (cseco ^ EcNRgammaseco))))))
  let lunggrowrateEcR (lunggrowthprop * (EcRxm * (1 - ((cprim ^ EcRgammaprim) / ((EcRec50prim ^ EcRgammaprim) + (cprim ^ EcRgammaprim)))) * (1 - ((cseco ^ EcRgammaseco) / ((EcRec50seco ^ EcRgammaseco) + (cseco ^ EcRgammaseco))))))
  let lunggrowrateSaNR (lunggrowthprop * (SaNRxm * (1 - ((cprim ^ SaNRgammaprim) / ((SaNRec50prim ^ SaNRgammaprim) + (cprim ^ SaNRgammaprim)))) * (1 - ((cseco ^ SaNRgammaseco) / ((SaNRec50seco ^ SaNRgammaseco) + (cseco ^ SaNRgammaseco))))))
  let lunggrowrateSaR (lunggrowthprop * (SaRxm * (1 - ((cprim ^ SaRgammaprim) / ((SaRec50prim ^ SaRgammaprim) + (cprim ^ SaRgammaprim)))) * (1 - ((cseco ^ SaRgammaseco) / ((SaRec50seco ^ SaRgammaseco) + (cseco ^ SaRgammaseco))))))

  ; for testing - grow lung pathogenic bacteria population according to growth rate only

  ;let lunggrowquantEcNR (lunggrowrateEcNR * lungquantEcNR)
  ;let lunggrowquantEcR (lunggrowrateEcR * lungquantEcR)
  ;let lunggrowquantSaNR (lunggrowrateSaNR * lungquantSaNR)
  ;let lunggrowquantSaR (lunggrowrateSaR * lungquantSaR)

  ; for testing - grow / diminish lung pathogenic bacteria population according to growth rate only innate immune system only

  ;let lunggrowquantEcNR (- ((lungmechclearrate * lungquantEcNR) / (lungphagoclearrate + lungquantEcNR)))
  ;let lunggrowquantEcR (- ((lungmechclearrate * lungquantEcR) / (lungphagoclearrate + lungquantEcR)))
  ;let lunggrowquantSaNR (- ((lungmechclearrate * lungquantSaNR) / (lungphagoclearrate + lungquantSaNR)))
  ;let lunggrowquantSaR (- ((lungmechclearrate * lungquantSaR) / (lungphagoclearrate + lungquantSaR)))

  ; for testing - grow using combination of growth rate and innate immune

  ;let lunggrowquantEcNR ((lunggrowrateEcNR * lungquantEcNR) - ((lungmechclearrate * lungquantEcNR) / (lungphagoclearrate + lungquantEcNR)))
  ;let lunggrowquantEcR ((lunggrowrateEcR * lungquantEcR) - ((lungmechclearrate * lungquantEcR) / (lungphagoclearrate + lungquantEcR)))
  ;let lunggrowquantSaNR ((lunggrowrateSaNR * lungquantSaNR) - ((lungmechclearrate * lungquantSaNR) / (lungphagoclearrate + lungquantSaNR)))
  ;let lunggrowquantSaR ((lunggrowrateSaR * lungquantSaR) - ((lungmechclearrate * lungquantSaR) / (lungphagoclearrate + lungquantSaR)))

  ; adjust phagorates by stage in life cycle

  let Ecphagorateadjusted Ecphagorate
  let Saphagorateadjusted Saphagorate

  if (stage = "starter") [
    set Ecphagorateadjusted (Ecphagorateadjusted * (stagephagorateadjust * 3))
    set Saphagorateadjusted (Saphagorateadjusted * (stagephagorateadjust * 3))
  ]

  if (stage = "finisher") [
    set Ecphagorateadjusted (Ecphagorateadjusted * stagephagorateadjust)
    set Saphagorateadjusted (Saphagorateadjusted * stagephagorateadjust)
  ]

  ; calculate growth of each lung pathogenic bacteria population - full model

  let lunggrowquantEcNR ((lunggrowrateEcNR * lungquantEcNR) - ((lungmechclearrate * lungquantEcNR) / (lungphagoclearrate + lungquantEcNR)) - ((Ecphagorateadjusted * lungquantImmune * lungquantEcNR) / (1 + (phagoinhib * lungquantEcNR))))
  let lunggrowquantEcR ((lunggrowrateEcR * lungquantEcR) - ((lungmechclearrate * lungquantEcR) / (lungphagoclearrate + lungquantEcR)) - ((Ecphagorateadjusted * lungquantImmune * lungquantEcR) / (1 + (phagoinhib * lungquantEcR))))
  let lunggrowquantSaNR ((lunggrowrateSaNR * lungquantSaNR) - ((lungmechclearrate * lungquantSaNR) / (lungphagoclearrate + lungquantSaNR)) - ((Saphagorateadjusted * lungquantImmune * lungquantSaNR) / (1 + (phagoinhib * lungquantSaNR))))
  let lunggrowquantSaR ((lunggrowrateSaR * lungquantSaR) - ((lungmechclearrate * lungquantSaR) / (lungphagoclearrate + lungquantSaR)) - ((Saphagorateadjusted * lungquantImmune * lungquantSaR) / (1 + (phagoinhib * lungquantSaR))))

  ;;;;;;;;;;;;;;;;;;; testing script

  ;if who = 0 [ask chicken 0 [print "prior quant" print lungquantEcNR]]

  ;if who = 0 [ask chicken 0 [print "result" print lunggrowquantEcNR]]

  ;let lunggrowterm1 (lunggrowrateEcNR * lungquantEcNR)
  ;if who = 0 [ask chicken 0 [print 1 print lunggrowterm1]]

  ;let lunggrowterm2 ((lungmechclearrate * lungquantEcNR) / (lungphagoclearrate + lungquantEcNR))
  ;if who = 0 [ask chicken 0 [print 2 print lunggrowterm2]]

  ;let lunggrowterm3 ((Ecphagorate * lungquantImmune * lungquantEcNR) / (1 + (phagoinhib * lungquantEcNR)))
  ;if who = 0 [ask chicken 0 [print 3 print lunggrowterm3]]

  ; modify each lung pathogenic bacteria population according to growth

  set lungquantEcNR (lungquantEcNR + lunggrowquantEcNR)
  set lungquantEcR (lungquantEcR + lunggrowquantEcR)
  set lungquantSaNR (lungquantSaNR + lunggrowquantSaNR)
  set lungquantSaR (lungquantSaR + lunggrowquantSaR)

  ; cutoff at 10 / limit at zero

  if lungquantEcNR < 10 [set lungquantEcNR 0]
  if lungquantEcR < 10 [set lungquantEcR 0]
  if lungquantSaNR < 10 [set lungquantSaNR 0]
  if lungquantSaR < 10 [set lungquantSaR 0]

  ;; set establish relevant lung population variables for use in damage and immune response

  let lungquantTotal (lungquantEcNR + lungquantEcR + lungquantSaNR + lungquantSaR)

  ;; add lung mech excretion to patches

  set patchEcNR (patchEcNR + ((lungmechclearrate * lungquantEcNR) / (lungphagoclearrate + lungquantEcNR)))
  set patchEcR (patchEcR + ((lungmechclearrate * lungquantEcR) / (lungphagoclearrate + lungquantEcR)))
  set patchSaNR (patchSaNR + ((lungmechclearrate * lungquantSaNR) / (lungphagoclearrate + lungquantSaNR)))
  set patchSaR (patchSaR + ((lungmechclearrate * lungquantSaR) / (lungphagoclearrate + lungquantSaR)))

  ;; calculate lung damage

  ; setup max 1 for increase coefficient term

  let damagecoeffmax 1

  let damagecoeff (damageincreaserate * lungquantTotal)

  if damagecoeff > 1 [set damagecoeff 1]

  ; apply model formula

  let damagechange ((damagecoeff * (1 - epiDamage)) - (damagerepairrate * epiDamage))

  ;if who = 0 [ask chicken 0 [print damagechange]]

  set epiDamage (epiDamage + damagechange)

  ;; calculate immune cells response, subject to infection threshold

  ; function below is literal from Mochan et al 2014 (Eq4)

  ;let lungimmunegrowquant ((- lungquantImmune + ((immuneinfluxrate * lungquantTotal) / (lungimmunesteady + lungquantTotal))) * (1 / immunedelay))

  ; function below applies coefficient to slow feedback loop between bacteria and immune response

    let lungimmunegrowquant 0.05 * ((- lungquantImmune + ((immuneinfluxrate * lungquantTotal) / (lungimmunesteady + lungquantTotal))) * (1 / immunedelay))

    ;if who = 0 [ask chicken 0 [print "signal" print lungimmunegrowquant]]

    ; save as immune growth signal for later use

    set immuneTracker fput lungimmunegrowquant but-last immuneTracker

    ;if who = 0 [ask chicken 0 [print immuneTracker]]

    ; update immune cells with prior signal

    set lungquantImmune (lungquantImmune + item (immunedelay - 1) immuneTracker)

    ;if who = 0 [ask chicken 0 [print "update" print lungquantImmune]]

  ; limit at steady state value

  if lungquantImmune < lungimmunesteady [set lungquantImmune lungimmunesteady]

end

; define procedure for strain growth on patches

to PatchStrainGrow

  let seasonadjust 1

  if season = "Continuous variable" [set seasonadjust seasonadjustobs]

  if season = "Early Dry" [set seasonadjust -20]
  if season = "Late Dry" [set seasonadjust 0.05]
  if season = "Early Rainy" [set seasonadjust 2]
  if season = "Late Rainy" [set seasonadjust 5]

  let patchgrowEcNR ((EcNRxm * patchEcNR) * patchgrowthprop * seasonadjust)
  let patchgrowEcR ((EcRxm * patchEcR) * patchgrowthprop * seasonadjust)
  let patchgrowSaNR ((SaNRxm * patchSaNR) * patchgrowthprop * seasonadjust)
  let patchgrowSaR ((SaRxm * patchSaR) * patchgrowthprop * seasonadjust)
  let patchgrowAuNR ((AuNRxm * patchSaNR) * patchgrowthprop * seasonadjust)
  let patchgrowAuR ((AuRxm * patchSaR) * patchgrowthprop * seasonadjust)

  set patchEcNR (patchEcNR + patchgrowEcNR)
  set patchEcR (patchEcR + patchgrowEcR)
  set patchSaNR (patchSaNR + patchgrowSaNR)
  set patchSaR (patchSaR + patchgrowSaR)
  set patchAuNR (patchAuNR + patchgrowAuNR)
  set patchAuR (patchAuR + patchgrowAuR)

end

; define procedure for update patch colors

to UpdateColors

  ifelse ((patchEcNR * 0.0015) >= 5.0E4) or ((patchEcR * 0.0015) >= 5.0E4) or ((patchSaNR * 0.0015) >= 5.0E4) or ((patchSaR * 0.0015) >= 5.0E4)

    [set pcolor 19 - (((patchEcNR + patchEcR + patchSaNR + patchSaR) / (1.0E9 * 0.5)) * 4)]

    [set pcolor 65 - (((patchEcNR + patchEcR + patchSaNR + patchSaR) / (1.0E9 * 0.5)) * 4)]

end

; define procedure for sigmoid weight gain cf. Mouffok et al 2019, nb sickness penalty / death -> 0

to GrowMass

  ; basic mass function

  let growX ((ticks / 24) / 7) - 3
  let gainVal (maxmass / (1 + e ^ (- growX)))

  let growX+1 (((ticks + 1) / 24) / 7) - 3
  let gainVal+1 (maxmass / (1 + e ^ (- growX+1)))

  let tickGain (gainVal+1 - gainVal)

  ; weight curve disregarding disease & death
  ;set mass (mass + tickGain)

  ; disease penalty - original, linked to lung damage

  ;let diseasepenalty (1 - epiDamage)

  ; disease penalty - new, linked to immune response

  let diseasepenalty 1 - ((lungquantImmune - lungimmunesteady) / (immunegrowthintersect - lungimmunesteady))

  ; add weight

  set mass (mass + (tickGain * diseasepenalty))

  ;; add death
  if viva = "Dead" [set mass 0]

end

; define procedure for treatment - different regimes

to TreatPrimHigh

  if AMUstance = "User choice" [

    ask chickens [set cprim 40 set countprim 72]

  ]

end

to TreatPrimMed

  if AMUstance = "User choice" [

    ask chickens [set cprim 10 set countprim 168]

  ]

end

to TreatPrimLow

  if AMUstance = "User choice" [

    ask chickens [set cprim 1 set countprim 336]

  ]

end

to TreatSeco

  if AMUstance = "User choice" [

    ask chickens [set cseco 40 set countseco 72]

  ]

end

; define procedures clean cage (nb now as adding fresh litter - effect is to bury existing bacteria, modelled as only 10% remain relevant to chickens [remainder still passed to disposal tracker for final outcomes])

to CleanCage

  if AMUstance = "User choice" [

    set clean-disposalPathoRsum (clean-disposalPathoRsum + (0.9 * sum [patchEcR] of patches) + (0.9 * sum [patchSaR] of patches))

    ask patches [set patchEcNR (0.1 * patchEcNR) set patchEcR (0.1 * patchEcR) set patchSaNR (0.1 * patchSaNR) set patchSaR (0.1 * patchSaR) set patchAuNR (0.1 * patchAuNR) set patchAuR (0.1 * patchAuR)]

    ask not-feeders [set pcolor 65]

    set clean-counter (clean-counter + 1)

  ]

end

to CleanCageAuto

    set clean-disposalPathoRsum (clean-disposalPathoRsum + (0.9 * sum [patchEcR] of patches) + (0.9 * sum [patchSaR] of patches))

    ask patches [set patchEcNR (0.1 * patchEcNR) set patchEcR (0.1 * patchEcR) set patchSaNR (0.1 * patchSaNR) set patchSaR (0.1 * patchSaR) set patchAuNR (0.1 * patchAuNR) set patchAuR (0.1 * patchAuR)]

    ask not-feeders [set pcolor 65]

    set clean-counter (clean-counter + 1)

end

; define procedure to dispose of dead chicken carcasses

to DisposeDead

  ifelse (lungquantEcNR + lungquantEcR) > (lungquantSaNR + lungquantSaR)
  [set deathcauseEc deathcauseEc + 1]
  [set deathcauseSa deathcauseSa + 1]

  set Deadcounter (Deadcounter + 1)

  carefully [set carcass-disposalPathoRsum (carcass-disposalPathoRsum + quantEcR + quantSaR + lungquantEcR + lungquantSaR)] [print error-message]

  die

end

; define procedure for movement & ChickenCount update

to MoveStep

  ; distinguish initial pen mobility

  if ticks <= 72
  [ask chickens [if viva = "Alive" [move-to one-of patch-set (neighbors with [member? self initialpen])]]]

  ; distribute chickens to feeders from initial pen

  if ticks = 73
  [ask chickens [if viva = "Alive" [move-to one-of [patches in-radius 6] of home-feeder]]]

  ; mobility over rest of life

  if ticks >= 74
  [ask chickens [if viva = "Alive" [

    ; stay close to home feeder (allowing chickens to swtich home feeder if closer)

    ifelse ticks mod 24 = 0
    [set home-feeder min-one-of feeders [distance myself] move-to one-of [patches in-radius 6] of home-feeder]

    ; lazy behaviour (move p set in interface by global variable moveprob)

    [if random-float 1 > (1 - moveprob) [move-to one-of neighbors]]
  ]]]

  ; update patch count of chickens

  ask patches [set chickenCount (count turtles-here)]

end

; define procedure to report outcomes

to Report-Outcomes

  ;output-print prim-tracker
  ;output-print seco-tracker
  ;output-print patchtotalsPathoR-tracker

  output-print (word "Strategy:" AMUstance)

  output-print ("")

  output-print (word "ECONOMIC OUTCOMES")

  output-print (word "1. Production")
  let meat sum [mass] of chickens
  output-print (word precision (meat / 1000) 2 " kg meat produced")
  output-print (word precision ((meat / ((count chickens + Deadcounter) * 2000)) * 100) 2 "% efficiency vs theoretical maximum")

  output-print ("")

  output-print (word "2. Income")
  let pricemeat 1.2
  output-print (word int (meat * pricemeat) " FCFA total income")
  output-print (word "(Per kilo price 1200 FCFA)")

  output-print ("")

  output-print (word "3. Expenses")
  let chicks (count chickens + Deadcounter) * 1000
  let feed ((sum [mass] of chickens + Deadmass) - ((count chickens + Deadcounter) * 200)) * 1.6 * 0.3
  let primcost sum prim-tracker * 0.01
  let secocost sum seco-tracker * 0.01
  let labour (clean-counter * 5000) + 10000
  output-print (word int (chicks + feed + primcost + secocost + labour) " FCFA total expenses")
  output-print (word "(Chick purchase plus feed over time plus medications plus labour)")
  output-print (word int chicks " on chicks at FCFA 1000 each")
  output-print (word int feed " on feed calculated as per gram of mass produced including birds that died and excluding start weight - FCR 1.6 and FCFA 300 per kilo (FCFA 15,000 per 50kg)")
  output-print (word "Medications costs included (assuming FCFA 0.01 to achieve 1 unit of concentration per bird in each case):")
  output-print (word int primcost " on Tetracycline")
  output-print (word int secocost " on fluoroquinolones")
  output-print (word int labour " on labour - fixed 10000 plus 5000 per adding of litter (" clean-counter " - times)")

  output-print ("")

  output-print (word "4. Profit")
  output-print (word int ((meat * pricemeat) - (chicks + feed + primcost + secocost + labour)) " FCFA gross profit")

  output-print ("")

  output-print (word "HUMAN HEALTH - RELATED OUTCOMES")

  output-print (word "1. AB residue in meat")
  let primresid (sum (sublist prim-tracker (ticks - withdrawmoment) (ticks)))
  let secoresid (sum (sublist seco-tracker (ticks - withdrawmoment) (ticks)))
  output-print (word int (primresid + secoresid) " units of antimicrobial medicine concentation in meat due to non-respect of withdrawal period of 10 days")

  output-print ("")

  output-print (word "2. Resistant pathogenic bacteria passed to meat supply chain")

  let finaltotalEcR (sum [quantEcR] of chickens + sum [lungquantEcR] of chickens)
  let finaltotalSaR (sum [quantSaR] of chickens + sum [lungquantSaR] of chickens)
  let finaltotalgutpathoR (sum [quantEcR] of chickens + sum [quantSaR] of chickens)
  let finaltotallungpathoR (sum [lungquantEcR] of chickens + sum [lungquantSaR] of chickens)

  output-print (word int (finaltotalEcR +  finaltotalSaR) " grand total")
  output-print (word "Including " int (finaltotalEcR) " resistant E Coli strain and " int (finaltotalSaR) " resistant Salmonella strain")
  output-print (word "Including " int (finaltotalgutpathoR) " in intestines and " int (finaltotallungpathoR) " in lungs")

  output-print ("")

  output-print (word "3. Exposure of producers & wider community to resistant pathogenic bacteria")
  output-print (word int (sum patchtotalsPathoR-tracker) " grand total patch EcR & SaR over whole production cycle")

  output-print ("")

  output-print (word "ENVIRONMENT - RELATED OUTCOMES")

  output-print (word "1. Waste disposal from cage cleaning - resistant pathogenic bacteria")
  output-print (word int clean-disposalPathoRsum " resistant pathogenic bacteria emitted to environment from final cleaning")

  output-print ("")

  output-print (word "2. Carcass disposal - resistant pathogenic bacteria")
  output-print (word int carcass-disposalPathoRsum " resistant pathogenic bacteria emitted to environment from carcass disposal")

end
@#$#@#$#@
GRAPHICS-WINDOW
223
104
1011
503
-1
-1
13.0
1
10
1
1
1
0
0
0
1
-30
29
-15
14
0
0
1
ticks
30.0

BUTTON
331
13
394
46
setup
setup
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
402
13
465
46
go
go
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

MONITOR
1113
186
1285
231
NIL
[quantEcNR] of chicken 0
17
1
11

MONITOR
1141
231
1285
276
NIL
[quantEcR] of chicken 0
17
1
11

MONITOR
1131
275
1285
320
NIL
[quantSaNR] of chicken 0
17
1
11

MONITOR
1139
320
1285
365
NIL
[quantSaR] of chicken 0
17
1
11

MONITOR
1130
365
1285
410
NIL
[quantAuNR] of chicken 0
17
1
11

MONITOR
1138
410
1285
455
NIL
[quantAuR] of chicken 0
17
1
11

PLOT
1084
34
1284
184
gut bacteria in chicken 0
NIL
NIL
0.0
10.0
0.0
10.0
true
false
"" ""
PENS
"default" 1.0 0 -16777216 true "" "plot [quantEcNR] of chicken 0"
"pen-1" 1.0 0 -7500403 true "" "plot [quantEcR] of chicken 0"
"pen-2" 1.0 0 -2674135 true "" "plot [quantSaNR] of chicken 0"
"pen-3" 1.0 0 -955883 true "" "plot [quantSaR] of chicken 0"
"pen-4" 1.0 0 -6459832 true "" "plot [quantAuNR] of chicken 0"
"pen-5" 1.0 0 -1184463 true "" "plot [quantAuR] of chicken 0"

PLOT
1077
539
1277
689
mass of chicken 0
NIL
NIL
0.0
10.0
0.0
10.0
true
false
"" ""
PENS
"default" 1.0 0 -16777216 true "" "plot [mass] of chicken 0"

MONITOR
527
10
584
55
Day
int (ticks / 24)
17
1
11

MONITOR
13
168
140
213
Example tetracycline
[cprim] of chicken 0
17
1
11

MONITOR
14
232
141
277
Example fluoroquinolone
[cseco] of chicken 0
17
1
11

BUTTON
12
86
130
119
Give high tetra.
TreatPrimHigh
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
18
113
131
146
Give medium tetra.
TreatPrimMed
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
25
140
130
173
Give low tetra.
TreatPrimLow
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
13
203
130
236
Give fluoroquinolone
TreatSeco\n
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
6
277
98
310
Add litter
CleanCage
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

MONITOR
961
645
1075
690
NIL
[viva] of chicken 0
17
1
11

PLOT
1284
34
1484
184
Lung bacteria in chicken 0
NIL
NIL
0.0
10.0
0.0
10.0
true
false
"" ""
PENS
"default" 1.0 0 -16777216 true "" "plot [lungquantEcNR] of chicken 0"
"pen-1" 1.0 0 -7500403 true "" "plot [lungquantEcR] of chicken 0"
"pen-2" 1.0 0 -2674135 true "" "plot [lungquantSaNR] of chicken 0"
"pen-3" 1.0 0 -955883 true "" "plot [lungquantSaR] of chicken 0"

PLOT
1286
380
1486
530
Lung immune cells in chicken 0
NIL
NIL
0.0
10.0
0.0
10.0
true
false
"" ""
PENS
"default" 1.0 0 -16777216 true "" "plot [lungquantImmune] of chicken 0\n"

MONITOR
6
311
172
356
NIL
[patchEcNR] of patch 10 10
17
1
11

MONITOR
6
356
164
401
NIL
[patchEcR] of patch 10 10
17
1
11

MONITOR
6
400
173
445
NIL
[patchSaNR] of patch 10 10
17
1
11

MONITOR
6
444
165
489
NIL
[patchSaR] of patch 10 10
17
1
11

PLOT
1286
539
1486
689
Lung damage in chicken 0
NIL
NIL
0.0
10.0
0.0
1.0
true
false
"" ""
PENS
"default" 1.0 0 -16777216 true "" "plot [epiDamage] of chicken 0"
"pen-1" 1.0 0 -2674135 true "" "plot 0.8"
"pen-2" 1.0 0 -5987164 true "" "plot 0.2"

MONITOR
1285
186
1460
231
NIL
[lungquantEcNR] of chicken 0
17
1
11

MONITOR
1285
231
1452
276
NIL
[lungquantEcR] of chicken 0
17
1
11

MONITOR
1285
276
1462
321
NIL
[lungquantSaNR] of chicken 0
17
1
11

MONITOR
1285
320
1454
365
NIL
[lungquantSaR] of chicken 0
17
1
11

MONITOR
369
557
446
602
NIL
Deadcounter
17
1
11

PLOT
463
546
828
696
total mass
NIL
NIL
0.0
10.0
0.0
10.0
true
true
"" ""
PENS
"actual" 1.0 0 -16777216 true "" "plot sum [mass] of chickens"
"potential" 1.0 0 -7500403 true "" "plot sum [maxmass] of chickens"

CHOOSER
6
495
164
540
season
season
"Continuous variable" "Early Dry" "Late Dry" "Early Rainy" "Late Rainy" "NA"
4

TEXTBOX
14
15
362
44
Implement simulation - setup & go
20
0.0
1

TEXTBOX
33
1403
183
1428
Set parameters
20
0.0
1

TEXTBOX
28
726
178
753
Outcomes
22
0.0
1

OUTPUT
129
766
1325
1386
11

INPUTBOX
47
1478
134
1538
EcNRxm_mean
0.18
1
0
Number

INPUTBOX
134
1478
201
1538
EcNRxm_sd
0.02
1
0
Number

INPUTBOX
47
1538
135
1598
EcRxm_mean
0.18
1
0
Number

INPUTBOX
134
1538
201
1598
EcRxm_sd
0.02
1
0
Number

INPUTBOX
200
1538
302
1598
EcRxm_maxfitcost%
30.0
1
0
Number

INPUTBOX
47
1606
136
1666
SaNRxm_mean
0.18
1
0
Number

INPUTBOX
136
1606
202
1666
SaNRxm_sd
0.02
1
0
Number

INPUTBOX
47
1666
137
1726
SaRxm_mean
0.18
1
0
Number

INPUTBOX
136
1666
202
1726
SaRxm_sd
0.02
1
0
Number

INPUTBOX
201
1666
302
1726
SaRxm_maxfitcost%
30.0
1
0
Number

INPUTBOX
48
1732
138
1792
AuNRxm_mean
0.18
1
0
Number

INPUTBOX
137
1732
203
1792
AuNRxm_sd
0.03
1
0
Number

INPUTBOX
48
1792
138
1852
AuRxm_mean
0.18
1
0
Number

INPUTBOX
137
1792
203
1852
AuRxm_sd
0.02
1
0
Number

INPUTBOX
202
1792
302
1852
AuRxm_maxfitcost%
30.0
1
0
Number

INPUTBOX
392
1478
467
1538
EcNRec50primmax
4.0
1
0
Number

INPUTBOX
318
1478
392
1538
EcNRec50primmin
0.1
1
0
Number

INPUTBOX
467
1478
545
1538
EcNRgammaprimmin
1.0
1
0
Number

INPUTBOX
544
1478
621
1538
EcNRgammaprimmax
3.0
1
0
Number

TEXTBOX
46
1451
196
1469
STRAIN GROWTH
11
0.0
1

TEXTBOX
316
1451
507
1469
STRAIN RESISTANCE - TETRACYCLINE
11
0.0
1

INPUTBOX
631
1478
708
1538
EcNRec50secomin
0.1
1
0
Number

INPUTBOX
708
1478
780
1538
EcNRec50secomax
4.0
1
0
Number

INPUTBOX
780
1478
855
1538
EcNRgammasecomin
1.0
1
0
Number

INPUTBOX
855
1478
929
1538
EcNRgammasecomax
3.0
1
0
Number

TEXTBOX
631
1453
830
1481
STRAIN RESISTANCE - QUINOLONES
11
0.0
1

TEXTBOX
48
1869
241
1897
PATCH EXCRETION & GROWTH
11
0.0
1

TEXTBOX
48
1870
198
1888
NIL
11
0.0
1

TEXTBOX
1012
1451
1162
1469
SEEDING
11
0.0
1

TEXTBOX
1210
1450
1360
1468
INTAKE
11
0.0
1

TEXTBOX
1355
1448
1505
1466
REMOVAL
11
0.0
1

TEXTBOX
286
1867
473
1895
INFECTION
11
0.0
1

TEXTBOX
675
1866
825
1884
IMMUNE SYSTEM
11
0.0
1

INPUTBOX
1354
1472
1452
1532
cutoffvalobs
5000000.0
1
0
Number

INPUTBOX
1354
1532
1453
1592
removprobobs
0.5
1
0
Number

INPUTBOX
51
1887
206
1947
excreterateobs
0.01
1
0
Number

INPUTBOX
1210
1473
1311
1533
intakerateobs
0.001
1
0
Number

TEXTBOX
1328
1867
1478
1885
LUNG DAMAGE
11
0.0
1

INPUTBOX
1010
1477
1146
1537
startcap
2.5E8
1
0
Number

INPUTBOX
1010
1540
1145
1600
propEc
0.2
1
0
Number

INPUTBOX
1010
1663
1148
1723
propSa
0.2
1
0
Number

INPUTBOX
1083
1599
1146
1659
maxRinEc
0.8
1
0
Number

INPUTBOX
1010
1600
1072
1660
minRinEc
0.5
1
0
Number

INPUTBOX
1010
1725
1073
1785
minRinSa
0.2
1
0
Number

INPUTBOX
1088
1725
1148
1785
maxRinSa
0.5
1
0
Number

INPUTBOX
1009
1788
1147
1848
RinAu
0.2
1
0
Number

INPUTBOX
51
1947
206
2007
patchgrowthpropobs
0.002
1
0
Number

INPUTBOX
287
1888
379
1948
infectprobsobs
0.01
1
0
Number

INPUTBOX
380
1888
482
1948
contamweightobs
1.0
1
0
Number

INPUTBOX
482
1888
587
1948
surroundweightobs
1.0
1
0
Number

INPUTBOX
287
1948
420
2008
infectscalesobs
6000000.0
1
0
Number

INPUTBOX
419
1948
503
2008
contamscaleobs
0.05
1
0
Number

INPUTBOX
501
1948
587
2008
surroundscaleobs
0.5
1
0
Number

INPUTBOX
452
2008
587
2068
surrounddelayobs
12.0
1
0
Number

INPUTBOX
673
1889
828
1949
lunggrowthpropobs
0.1
1
0
Number

INPUTBOX
673
1949
828
2009
lungmechclearrateobs
1000.0
1
0
Number

INPUTBOX
673
2008
828
2068
lungphagoclearrateobs
1000000.0
1
0
Number

INPUTBOX
828
1889
983
1949
lungimmunesteadyobs
1000.0
1
0
Number

INPUTBOX
828
1949
907
2009
Ecphagorateobs
1.0E-7
1
0
Number

INPUTBOX
907
1949
983
2009
Saphagorateobs
5.0E-5
1
0
Number

INPUTBOX
983
1889
1138
1949
phagoinhibobs
1.0E-13
1
0
Number

INPUTBOX
983
1949
1138
2009
immuneinfluxrateobs
10000.0
1
0
Number

INPUTBOX
983
2008
1138
2068
immunedelayobs
6.0
1
0
Number

INPUTBOX
1326
1889
1481
1949
damageincreaserateobs
1.0E-11
1
0
Number

INPUTBOX
1326
1949
1481
2009
damagerepairrateobs
0.1
1
0
Number

MONITOR
600
10
678
55
Stage
stage
17
1
11

INPUTBOX
828
2009
983
2069
stagephagorateadjustobs
0.7
1
0
Number

PLOT
23
545
351
695
Apparent infections in live flock
NIL
%
0.0
10.0
0.0
1.0
true
true
"" ""
PENS
"Colibaccillosis" 1.0 0 -2674135 true "" "plot (apparentEc / (count chickens) * 100)"
"Salmonellosis" 1.0 0 -14070903 true "" "plot (apparentSa / count chickens * 100)"

MONITOR
369
602
426
647
%Ec
precision (deathcauseEc / Deadcounter * 100) 2
17
1
11

MONITOR
369
646
426
691
%Sa
precision (deathcauseSa / Deadcounter * 100) 2
17
1
11

CHOOSER
12
40
252
85
AMUstance
AMUstance
"User choice" "No intervention" "Responsive quick stop" "Responsive complete" "Prophylactic" "Clean responsive" "Regular clean" "Responsive quick stop + clean" "Responsive complete + reg. clean"
0

INPUTBOX
318
1537
392
1597
EcRec50primmin
16.0
1
0
Number

INPUTBOX
392
1537
468
1597
EcRec50primmax
200.0
1
0
Number

INPUTBOX
467
1537
545
1597
EcRgammaprimmin
8.0
1
0
Number

INPUTBOX
544
1537
621
1597
EcRgammaprimmax
20.0
1
0
Number

INPUTBOX
631
1538
709
1598
EcRec50secomin
8.0
1
0
Number

INPUTBOX
708
1538
780
1598
EcRec50secomax
32.0
1
0
Number

INPUTBOX
780
1538
856
1598
EcRgammasecomin
2.0
1
0
Number

INPUTBOX
856
1538
929
1598
EcRgammasecomax
5.0
1
0
Number

INPUTBOX
318
1605
393
1665
SaNRec50primmin
0.1
1
0
Number

INPUTBOX
393
1605
468
1665
SaNRec50primmax
4.0
1
0
Number

INPUTBOX
468
1605
546
1665
SaNRgammaprimmin
1.0
1
0
Number

INPUTBOX
546
1605
621
1665
SaNRgammaprimmax
3.0
1
0
Number

INPUTBOX
631
1606
708
1666
SaNRec50secomin
0.1
1
0
Number

INPUTBOX
708
1606
781
1666
SaNRec50secomax
4.0
1
0
Number

INPUTBOX
781
1606
856
1666
SaNRgammasecomin
1.0
1
0
Number

INPUTBOX
856
1606
929
1666
SaNRgammasecomax
3.0
1
0
Number

INPUTBOX
318
1665
393
1725
SaRec50primmin
8.0
1
0
Number

INPUTBOX
392
1665
469
1725
SaRec50primmax
100.0
1
0
Number

INPUTBOX
469
1665
547
1725
SaRgammaprimmin
4.0
1
0
Number

INPUTBOX
547
1665
621
1725
SaRgammaprimmax
10.0
1
0
Number

INPUTBOX
631
1665
708
1725
SaRec50secomin
1.0
1
0
Number

INPUTBOX
708
1665
782
1725
SaRec50secomax
128.0
1
0
Number

INPUTBOX
782
1666
857
1726
SaRgammasecomin
4.0
1
0
Number

INPUTBOX
856
1666
929
1726
SaRgammasecomax
10.0
1
0
Number

INPUTBOX
318
1733
393
1793
AuNRec50primmin
0.1
1
0
Number

INPUTBOX
393
1733
469
1793
AuNRec50primmax
4.0
1
0
Number

INPUTBOX
469
1734
548
1794
AuNRgammaprimmin
1.0
1
0
Number

INPUTBOX
547
1734
621
1794
AuNRgammaprimmax
3.0
1
0
Number

INPUTBOX
631
1734
708
1794
AuNRec50secomin
0.1
1
0
Number

INPUTBOX
707
1734
782
1794
AuNRec50secomax
4.0
1
0
Number

INPUTBOX
781
1734
857
1794
AuNRgammasecomin
1.0
1
0
Number

INPUTBOX
856
1734
929
1794
AuNRgammasecomax
3.0
1
0
Number

INPUTBOX
318
1793
393
1853
AuRec50primmin
16.0
1
0
Number

INPUTBOX
393
1793
470
1853
AuRec50primmax
200.0
1
0
Number

INPUTBOX
469
1793
548
1853
AuRgammaprimmin
8.0
1
0
Number

INPUTBOX
547
1793
621
1853
AuRgammaprimmax
20.0
1
0
Number

INPUTBOX
631
1793
708
1853
AuRec50secomin
4.0
1
0
Number

INPUTBOX
707
1793
781
1853
AuRec50secomax
50.0
1
0
Number

INPUTBOX
781
1793
856
1853
AuRgammasecomin
2.0
1
0
Number

INPUTBOX
854
1793
929
1853
AuRgammasecomax
5.0
1
0
Number

MONITOR
912
597
1075
642
NIL
[home-feeder] of chicken 0
17
1
11

MONITOR
956
549
1075
594
NIL
[mass] of chicken 0
17
1
11

MONITOR
778
650
828
695
% live chickens at slaughter weight
int ((count chickens with [mass > maxmass]) / (count chickens) * 100)
17
1
11

INPUTBOX
1151
1890
1306
1950
immunegrowthintersectobs
50000.0
1
0
Number

TEXTBOX
1155
1866
1305
1884
GROWTH
11
0.0
1

INPUTBOX
51
2006
206
2066
seasonadjustobs
-20.0
1
0
Number

TEXTBOX
1214
1655
1364
1673
MOBILITY
11
0.0
1

INPUTBOX
1212
1675
1367
1735
moveprob
0.3
1
0
Number

MONITOR
1328
766
1498
811
Mortality calibration
mortcalib
17
1
11

@#$#@#$#@
## WHAT IS IT?

(a general understanding of what the model is trying to show or explain)

## HOW IT WORKS

(what rules the agents use to create the overall behavior of the model)

## HOW TO USE IT

(how to use the model, including a description of each of the items in the Interface tab)

## THINGS TO NOTICE

(suggested things for the user to notice while running the model)

## THINGS TO TRY

(suggested things for the user to try to do (move sliders, switches, etc.) with the model)

## EXTENDING THE MODEL

(suggested things to add or change in the Code tab to make the model more complicated, detailed, accurate, etc.)

## NETLOGO FEATURES

(interesting or unusual features of NetLogo that the model uses, particularly in the Code tab; or where workarounds were needed for missing features)

## RELATED MODELS

(models in the NetLogo Models Library and elsewhere which are of related interest)

## CREDITS AND REFERENCES

(a reference to the model's URL on the web if it has one, as well as any other necessary credits, citations, and links)
@#$#@#$#@
default
true
0
Polygon -7500403 true true 150 5 40 250 150 205 260 250

airplane
true
0
Polygon -7500403 true true 150 0 135 15 120 60 120 105 15 165 15 195 120 180 135 240 105 270 120 285 150 270 180 285 210 270 165 240 180 180 285 195 285 165 180 105 180 60 165 15

arrow
true
0
Polygon -7500403 true true 150 0 0 150 105 150 105 293 195 293 195 150 300 150

box
false
0
Polygon -7500403 true true 150 285 285 225 285 75 150 135
Polygon -7500403 true true 150 135 15 75 150 15 285 75
Polygon -7500403 true true 15 75 15 225 150 285 150 135
Line -16777216 false 150 285 150 135
Line -16777216 false 150 135 15 75
Line -16777216 false 150 135 285 75

bug
true
0
Circle -7500403 true true 96 182 108
Circle -7500403 true true 110 127 80
Circle -7500403 true true 110 75 80
Line -7500403 true 150 100 80 30
Line -7500403 true 150 100 220 30

butterfly
true
0
Polygon -7500403 true true 150 165 209 199 225 225 225 255 195 270 165 255 150 240
Polygon -7500403 true true 150 165 89 198 75 225 75 255 105 270 135 255 150 240
Polygon -7500403 true true 139 148 100 105 55 90 25 90 10 105 10 135 25 180 40 195 85 194 139 163
Polygon -7500403 true true 162 150 200 105 245 90 275 90 290 105 290 135 275 180 260 195 215 195 162 165
Polygon -16777216 true false 150 255 135 225 120 150 135 120 150 105 165 120 180 150 165 225
Circle -16777216 true false 135 90 30
Line -16777216 false 150 105 195 60
Line -16777216 false 150 105 105 60

car
false
0
Polygon -7500403 true true 300 180 279 164 261 144 240 135 226 132 213 106 203 84 185 63 159 50 135 50 75 60 0 150 0 165 0 225 300 225 300 180
Circle -16777216 true false 180 180 90
Circle -16777216 true false 30 180 90
Polygon -16777216 true false 162 80 132 78 134 135 209 135 194 105 189 96 180 89
Circle -7500403 true true 47 195 58
Circle -7500403 true true 195 195 58

circle
false
0
Circle -7500403 true true 0 0 300

circle 2
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240

cow
false
0
Polygon -7500403 true true 200 193 197 249 179 249 177 196 166 187 140 189 93 191 78 179 72 211 49 209 48 181 37 149 25 120 25 89 45 72 103 84 179 75 198 76 252 64 272 81 293 103 285 121 255 121 242 118 224 167
Polygon -7500403 true true 73 210 86 251 62 249 48 208
Polygon -7500403 true true 25 114 16 195 9 204 23 213 25 200 39 123

cylinder
false
0
Circle -7500403 true true 0 0 300

dot
false
0
Circle -7500403 true true 90 90 120

face happy
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 255 90 239 62 213 47 191 67 179 90 203 109 218 150 225 192 218 210 203 227 181 251 194 236 217 212 240

face neutral
false
0
Circle -7500403 true true 8 7 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Rectangle -16777216 true false 60 195 240 225

face sad
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 168 90 184 62 210 47 232 67 244 90 220 109 205 150 198 192 205 210 220 227 242 251 229 236 206 212 183

fish
false
0
Polygon -1 true false 44 131 21 87 15 86 0 120 15 150 0 180 13 214 20 212 45 166
Polygon -1 true false 135 195 119 235 95 218 76 210 46 204 60 165
Polygon -1 true false 75 45 83 77 71 103 86 114 166 78 135 60
Polygon -7500403 true true 30 136 151 77 226 81 280 119 292 146 292 160 287 170 270 195 195 210 151 212 30 166
Circle -16777216 true false 215 106 30

flag
false
0
Rectangle -7500403 true true 60 15 75 300
Polygon -7500403 true true 90 150 270 90 90 30
Line -7500403 true 75 135 90 135
Line -7500403 true 75 45 90 45

flower
false
0
Polygon -10899396 true false 135 120 165 165 180 210 180 240 150 300 165 300 195 240 195 195 165 135
Circle -7500403 true true 85 132 38
Circle -7500403 true true 130 147 38
Circle -7500403 true true 192 85 38
Circle -7500403 true true 85 40 38
Circle -7500403 true true 177 40 38
Circle -7500403 true true 177 132 38
Circle -7500403 true true 70 85 38
Circle -7500403 true true 130 25 38
Circle -7500403 true true 96 51 108
Circle -16777216 true false 113 68 74
Polygon -10899396 true false 189 233 219 188 249 173 279 188 234 218
Polygon -10899396 true false 180 255 150 210 105 210 75 240 135 240

house
false
0
Rectangle -7500403 true true 45 120 255 285
Rectangle -16777216 true false 120 210 180 285
Polygon -7500403 true true 15 120 150 15 285 120
Line -16777216 false 30 120 270 120

leaf
false
0
Polygon -7500403 true true 150 210 135 195 120 210 60 210 30 195 60 180 60 165 15 135 30 120 15 105 40 104 45 90 60 90 90 105 105 120 120 120 105 60 120 60 135 30 150 15 165 30 180 60 195 60 180 120 195 120 210 105 240 90 255 90 263 104 285 105 270 120 285 135 240 165 240 180 270 195 240 210 180 210 165 195
Polygon -7500403 true true 135 195 135 240 120 255 105 255 105 285 135 285 165 240 165 195

line
true
0
Line -7500403 true 150 0 150 300

line half
true
0
Line -7500403 true 150 0 150 150

pentagon
false
0
Polygon -7500403 true true 150 15 15 120 60 285 240 285 285 120

person
false
0
Circle -7500403 true true 110 5 80
Polygon -7500403 true true 105 90 120 195 90 285 105 300 135 300 150 225 165 300 195 300 210 285 180 195 195 90
Rectangle -7500403 true true 127 79 172 94
Polygon -7500403 true true 195 90 240 150 225 180 165 105
Polygon -7500403 true true 105 90 60 150 75 180 135 105

plant
false
0
Rectangle -7500403 true true 135 90 165 300
Polygon -7500403 true true 135 255 90 210 45 195 75 255 135 285
Polygon -7500403 true true 165 255 210 210 255 195 225 255 165 285
Polygon -7500403 true true 135 180 90 135 45 120 75 180 135 210
Polygon -7500403 true true 165 180 165 210 225 180 255 120 210 135
Polygon -7500403 true true 135 105 90 60 45 45 75 105 135 135
Polygon -7500403 true true 165 105 165 135 225 105 255 45 210 60
Polygon -7500403 true true 135 90 120 45 150 15 180 45 165 90

sheep
false
15
Circle -1 true true 203 65 88
Circle -1 true true 70 65 162
Circle -1 true true 150 105 120
Polygon -7500403 true false 218 120 240 165 255 165 278 120
Circle -7500403 true false 214 72 67
Rectangle -1 true true 164 223 179 298
Polygon -1 true true 45 285 30 285 30 240 15 195 45 210
Circle -1 true true 3 83 150
Rectangle -1 true true 65 221 80 296
Polygon -1 true true 195 285 210 285 210 240 240 210 195 210
Polygon -7500403 true false 276 85 285 105 302 99 294 83
Polygon -7500403 true false 219 85 210 105 193 99 201 83

square
false
0
Rectangle -7500403 true true 30 30 270 270

square 2
false
0
Rectangle -7500403 true true 30 30 270 270
Rectangle -16777216 true false 60 60 240 240

star
false
0
Polygon -7500403 true true 151 1 185 108 298 108 207 175 242 282 151 216 59 282 94 175 3 108 116 108

target
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240
Circle -7500403 true true 60 60 180
Circle -16777216 true false 90 90 120
Circle -7500403 true true 120 120 60

tree
false
0
Circle -7500403 true true 118 3 94
Rectangle -6459832 true false 120 195 180 300
Circle -7500403 true true 65 21 108
Circle -7500403 true true 116 41 127
Circle -7500403 true true 45 90 120
Circle -7500403 true true 104 74 152

triangle
false
0
Polygon -7500403 true true 150 30 15 255 285 255

triangle 2
false
0
Polygon -7500403 true true 150 30 15 255 285 255
Polygon -16777216 true false 151 99 225 223 75 224

truck
false
0
Rectangle -7500403 true true 4 45 195 187
Polygon -7500403 true true 296 193 296 150 259 134 244 104 208 104 207 194
Rectangle -1 true false 195 60 195 105
Polygon -16777216 true false 238 112 252 141 219 141 218 112
Circle -16777216 true false 234 174 42
Rectangle -7500403 true true 181 185 214 194
Circle -16777216 true false 144 174 42
Circle -16777216 true false 24 174 42
Circle -7500403 false true 24 174 42
Circle -7500403 false true 144 174 42
Circle -7500403 false true 234 174 42

turtle
true
0
Polygon -10899396 true false 215 204 240 233 246 254 228 266 215 252 193 210
Polygon -10899396 true false 195 90 225 75 245 75 260 89 269 108 261 124 240 105 225 105 210 105
Polygon -10899396 true false 105 90 75 75 55 75 40 89 31 108 39 124 60 105 75 105 90 105
Polygon -10899396 true false 132 85 134 64 107 51 108 17 150 2 192 18 192 52 169 65 172 87
Polygon -10899396 true false 85 204 60 233 54 254 72 266 85 252 107 210
Polygon -7500403 true true 119 75 179 75 209 101 224 135 220 225 175 261 128 261 81 224 74 135 88 99

wheel
false
0
Circle -7500403 true true 3 3 294
Circle -16777216 true false 30 30 240
Line -7500403 true 150 285 150 15
Line -7500403 true 15 150 285 150
Circle -7500403 true true 120 120 60
Line -7500403 true 216 40 79 269
Line -7500403 true 40 84 269 221
Line -7500403 true 40 216 269 79
Line -7500403 true 84 40 221 269

wolf
false
0
Polygon -16777216 true false 253 133 245 131 245 133
Polygon -7500403 true true 2 194 13 197 30 191 38 193 38 205 20 226 20 257 27 265 38 266 40 260 31 253 31 230 60 206 68 198 75 209 66 228 65 243 82 261 84 268 100 267 103 261 77 239 79 231 100 207 98 196 119 201 143 202 160 195 166 210 172 213 173 238 167 251 160 248 154 265 169 264 178 247 186 240 198 260 200 271 217 271 219 262 207 258 195 230 192 198 210 184 227 164 242 144 259 145 284 151 277 141 293 140 299 134 297 127 273 119 270 105
Polygon -7500403 true true -1 195 14 180 36 166 40 153 53 140 82 131 134 133 159 126 188 115 227 108 236 102 238 98 268 86 269 92 281 87 269 103 269 113

x
false
0
Polygon -7500403 true true 270 75 225 30 30 225 75 270
Polygon -7500403 true true 30 75 75 30 270 225 225 270
@#$#@#$#@
NetLogo 6.2.2
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
<experiments>
  <experiment name="parameter space exploration - infection &amp; immunity" repetitions="4" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="2000"/>
    <metric>count chickens</metric>
    <metric>sum [mass] of chickens</metric>
    <enumeratedValueSet variable="infectprobsobs">
      <value value="0.001"/>
      <value value="0.005"/>
      <value value="0.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="contamweightobs">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="surroundweightobs">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="infectscalesobs">
      <value value="1000000"/>
      <value value="4000000"/>
      <value value="10000000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="contamscaleobs">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="surroundscaleobs">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="lunggrowthpropobs">
      <value value="0.1"/>
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="lungmechclearrateobs">
      <value value="5000"/>
      <value value="10000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="lungimmunesteadyobs">
      <value value="1000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="immuneinfluxrateobs">
      <value value="50000"/>
      <value value="100000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="lungphagoclearrateobs">
      <value value="5000000"/>
      <value value="10000000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Ecphagorateobs">
      <value value="1.0E-7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Saphagorateobs">
      <value value="5.0E-4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stagephagorateadjustobs">
      <value value="0.7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="phagoinhibobs">
      <value value="1.0E-11"/>
      <value value="1.0E-10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="damageincreaserateobs">
      <value value="1.0E-10"/>
      <value value="1.0E-11"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="damagerepairrateobs">
      <value value="0.04"/>
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="immunegrowthintersectobs">
      <value value="50000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="season">
      <value value="&quot;Late Rainy&quot;"/>
      <value value="&quot;Early Dry&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="immunedelayobs">
      <value value="6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="surrounddelayobs">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AuRec50primmin">
      <value value="16"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="SaNRec50secomax">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="SaNRxm_mean">
      <value value="0.18"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AuNRgammasecomax">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="SaRgammasecomin">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="EcNRgammaprimmax">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="SaNRgammasecomax">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="EcNRgammasecomin">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AuRgammaprimmax">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AuNRec50secomin">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="EcRec50primmax">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="maxRinEc">
      <value value="0.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="EcRxm_mean">
      <value value="0.18"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AuNRxm_sd">
      <value value="0.03"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="EcNRec50secomax">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="SaRgammaprimmax">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="excreterateobs">
      <value value="0.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AuRec50secomax">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="minRinSa">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="SaNRec50secomin">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AuNRgammasecomin">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AMUstance">
      <value value="&quot;User choice&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="SaRec50primmax">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AuNRec50primmax">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AuRgammaprimmin">
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="SaRxm_mean">
      <value value="0.18"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="SaRec50secomin">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="EcRec50primmin">
      <value value="16"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AuRxm_sd">
      <value value="0.02"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="EcRgammaprimmax">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="SaNRec50primmax">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="EcRgammasecomin">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AuNRgammaprimmax">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="SaRgammaprimmin">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="intakerateobs">
      <value value="0.001"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="SaNRgammaprimmax">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="EcNRgammaprimmin">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="SaNRgammasecomin">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cutoffvalobs">
      <value value="5000000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="EcRxm_maxfitcost%">
      <value value="30"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AuNRec50primmin">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="propSa">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="EcNRec50primmax">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="SaRxm_maxfitcost%">
      <value value="30"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="SaRxm_sd">
      <value value="0.02"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AuRec50primmax">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="SaNRxm_sd">
      <value value="0.02"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="EcNRec50secomin">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="SaNRec50primmin">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AuNRgammaprimmin">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AuRec50secomin">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="EcNRxm_sd">
      <value value="0.02"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="EcNRgammasecomax">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="EcNRxm_mean">
      <value value="0.18"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="patchgrowthpropobs">
      <value value="0.002"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="removprobobs">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AuRgammasecomax">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="EcRxm_sd">
      <value value="0.02"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="SaRec50primmin">
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="maxRinSa">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="EcRec50secomax">
      <value value="32"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="EcRgammaprimmin">
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="SaRgammasecomax">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AuRxm_maxfitcost%">
      <value value="30"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="minRinEc">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="SaNRgammaprimmin">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AuNRec50secomax">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="SaRec50secomax">
      <value value="128"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="startcap">
      <value value="250000000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AuNRxm_mean">
      <value value="0.18"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AuRgammasecomin">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="propEc">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AuRxm_mean">
      <value value="0.18"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="EcRec50secomin">
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="EcRgammasecomax">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="EcNRec50primmin">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="RinAu">
      <value value="0.2"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="AMUstances-experiments" repetitions="100" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="2000"/>
    <metric>count chickens</metric>
    <metric>sum [mass] of chickens</metric>
    <enumeratedValueSet variable="AMUstance">
      <value value="&quot;No intervention&quot;"/>
      <value value="&quot;Responsive quick stop&quot;"/>
      <value value="&quot;Responsive complete&quot;"/>
      <value value="&quot;Prophylactic&quot;"/>
      <value value="&quot;Clean responsive&quot;"/>
      <value value="&quot;Regular clean&quot;"/>
      <value value="&quot;Responsive quick stop + clean&quot;"/>
      <value value="&quot;Responsive complete + reg. clean&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="season">
      <value value="&quot;Continuous variable&quot;"/>
      <value value="&quot;Early Dry&quot;"/>
      <value value="&quot;Late Dry&quot;"/>
      <value value="&quot;Early Rainy&quot;"/>
      <value value="&quot;Late Rainy&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="surroundscaleobs">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="lungphagoclearrateobs">
      <value value="1000000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AuRec50primmin">
      <value value="16"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="SaNRec50secomax">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="SaNRxm_mean">
      <value value="0.18"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AuNRgammasecomax">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="SaRgammasecomin">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="EcNRgammaprimmax">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="SaNRgammasecomax">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="EcNRgammasecomin">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="damagerepairrateobs">
      <value value="0.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AuRgammaprimmax">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="immunegrowthintersectobs">
      <value value="50000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AuNRec50secomin">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="EcRec50primmax">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="maxRinEc">
      <value value="0.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="lungimmunesteadyobs">
      <value value="1000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="EcRxm_mean">
      <value value="0.18"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AuNRxm_sd">
      <value value="0.03"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="EcNRec50secomax">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="SaRgammaprimmax">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AuRec50secomax">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="excreterateobs">
      <value value="0.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="immuneinfluxrateobs">
      <value value="10000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="SaNRec50secomin">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="minRinSa">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AuNRgammasecomin">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="SaRec50primmax">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AuNRec50primmax">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AuRgammaprimmin">
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="SaRxm_mean">
      <value value="0.18"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="SaRec50secomin">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="EcRec50primmin">
      <value value="16"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AuRxm_sd">
      <value value="0.02"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="EcRgammaprimmax">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="SaNRec50primmax">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="EcRgammasecomin">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AuNRgammaprimmax">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="SaRgammaprimmin">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="intakerateobs">
      <value value="0.001"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="SaNRgammaprimmax">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="EcNRgammaprimmin">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="SaNRgammasecomin">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cutoffvalobs">
      <value value="5000000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Ecphagorateobs">
      <value value="1.0E-7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="EcRxm_maxfitcost%">
      <value value="30"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AuNRec50primmin">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="damageincreaserateobs">
      <value value="5.0E-9"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="propSa">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="EcNRec50primmax">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="SaRxm_sd">
      <value value="0.02"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="SaRxm_maxfitcost%">
      <value value="30"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stagephagorateadjustobs">
      <value value="0.7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AuRec50primmax">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="SaNRxm_sd">
      <value value="0.02"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="EcNRec50secomin">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="SaNRec50primmin">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AuNRgammaprimmin">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AuRec50secomin">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="EcNRxm_sd">
      <value value="0.02"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="lunggrowthpropobs">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="EcNRgammasecomax">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="infectprobsobs">
      <value value="0.008"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="EcNRxm_mean">
      <value value="0.18"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AuRgammasecomax">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="removprobobs">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="EcRxm_sd">
      <value value="0.02"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="patchgrowthpropobs">
      <value value="0.0012"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="contamscaleobs">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="SaRec50primmin">
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="maxRinSa">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="EcRec50secomax">
      <value value="32"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="infectscalesobs">
      <value value="400000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="EcRgammaprimmin">
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="immunedelayobs">
      <value value="6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AuRxm_maxfitcost%">
      <value value="30"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="SaRgammasecomax">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="phagoinhibobs">
      <value value="1.0E-13"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="surrounddelayobs">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="minRinEc">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Saphagorateobs">
      <value value="5.0E-5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="SaNRgammaprimmin">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="contamweightobs">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="lungmechclearrateobs">
      <value value="1000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="SaRec50secomax">
      <value value="128"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AuNRec50secomax">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="startcap">
      <value value="250000000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AuNRxm_mean">
      <value value="0.18"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AuRgammasecomin">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="propEc">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AuRxm_mean">
      <value value="0.18"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="surroundweightobs">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="EcRec50secomin">
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="seasonadjustobs">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="moveprob">
      <value value="0.25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="EcNRec50primmin">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="EcRgammasecomax">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="RinAu">
      <value value="0.2"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="nointerv-fitting" repetitions="20" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="2000"/>
    <metric>count(chickens)</metric>
    <metric>deadcounter</metric>
    <enumeratedValueSet variable="season">
      <value value="&quot;Continuous variable&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="seasonadjustobs">
      <value value="0"/>
      <value value="-10"/>
      <value value="-20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="patchgrowthpropobs">
      <value value="0.0012"/>
      <value value="0.002"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="lunggrowthpropobs">
      <value value="0.09"/>
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="infectscalesobs">
      <value value="100000"/>
      <value value="300000"/>
      <value value="500000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="damageincreaserateobs">
      <value value="3.0E-9"/>
      <value value="5.0E-9"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="damagerepairrateobs">
      <value value="0.01"/>
      <value value="0.02"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="moveprob">
      <value value="0.1"/>
      <value value="0.3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="surroundscaleobs">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="lungphagoclearrateobs">
      <value value="1000000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AuRec50primmin">
      <value value="16"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="SaNRec50secomax">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="SaNRxm_mean">
      <value value="0.18"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AuNRgammasecomax">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="EcNRgammaprimmax">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="SaRgammasecomin">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="SaNRgammasecomax">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="EcNRgammasecomin">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AuRgammaprimmax">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="immunegrowthintersectobs">
      <value value="50000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AuNRec50secomin">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="EcRec50primmax">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="maxRinEc">
      <value value="0.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="lungimmunesteadyobs">
      <value value="1000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="EcRxm_mean">
      <value value="0.18"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AuNRxm_sd">
      <value value="0.03"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="EcNRec50secomax">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="SaRgammaprimmax">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="excreterateobs">
      <value value="0.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AuRec50secomax">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="immuneinfluxrateobs">
      <value value="10000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="SaNRec50secomin">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="minRinSa">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AuNRgammasecomin">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AMUstance">
      <value value="&quot;User choice&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AuNRec50primmax">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="SaRec50primmax">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AuRgammaprimmin">
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="SaRxm_mean">
      <value value="0.18"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="SaRec50secomin">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="EcRec50primmin">
      <value value="16"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AuRxm_sd">
      <value value="0.02"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="EcRgammaprimmax">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="SaNRec50primmax">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="EcRgammasecomin">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AuNRgammaprimmax">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="SaRgammaprimmin">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="intakerateobs">
      <value value="0.001"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="SaNRgammaprimmax">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="EcNRgammaprimmin">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="SaNRgammasecomin">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cutoffvalobs">
      <value value="5000000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Ecphagorateobs">
      <value value="1.0E-7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="EcRxm_maxfitcost%">
      <value value="30"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AuNRec50primmin">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="damageincreaserateobs">
      <value value="5.0E-9"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="propSa">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="EcNRec50primmax">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="SaRxm_sd">
      <value value="0.02"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="SaRxm_maxfitcost%">
      <value value="30"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stagephagorateadjustobs">
      <value value="0.7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AuRec50primmax">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="SaNRxm_sd">
      <value value="0.02"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="EcNRec50secomin">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="SaNRec50primmin">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AuNRgammaprimmin">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AuRec50secomin">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="EcNRxm_sd">
      <value value="0.02"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="EcNRgammasecomax">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="infectprobsobs">
      <value value="0.008"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="EcNRxm_mean">
      <value value="0.18"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AuRgammasecomax">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="removprobobs">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="EcRxm_sd">
      <value value="0.02"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="contamscaleobs">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="SaRec50primmin">
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="maxRinSa">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="EcRec50secomax">
      <value value="32"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="EcRgammaprimmin">
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="infectscalesobs">
      <value value="400000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="immunedelayobs">
      <value value="6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="SaRgammasecomax">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AuRxm_maxfitcost%">
      <value value="30"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="phagoinhibobs">
      <value value="1.0E-13"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="surrounddelayobs">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="minRinEc">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Saphagorateobs">
      <value value="5.0E-5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="SaNRgammaprimmin">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="contamweightobs">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="lungmechclearrateobs">
      <value value="1000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AuNRec50secomax">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="SaRec50secomax">
      <value value="128"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="startcap">
      <value value="250000000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AuNRxm_mean">
      <value value="0.18"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AuRgammasecomin">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="propEc">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AuRxm_mean">
      <value value="0.18"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="EcRec50secomin">
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="surroundweightobs">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="EcNRec50primmin">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="EcRgammasecomax">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="RinAu">
      <value value="0.2"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="calibration13.4" repetitions="20" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="2000"/>
    <metric>deadcounter</metric>
    <enumeratedValueSet variable="seasonadjustobs">
      <value value="-20"/>
      <value value="-10"/>
      <value value="0"/>
      <value value="10"/>
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="infectprobsobs">
      <value value="0.01"/>
      <value value="0.011"/>
      <value value="0.012"/>
      <value value="0.013"/>
      <value value="0.014"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="infectscalesobs">
      <value value="9000000"/>
      <value value="10000000"/>
      <value value="11000000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="lunggrowthpropobs">
      <value value="0.1"/>
      <value value="0.11"/>
      <value value="0.12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="surroundscaleobs">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="lungphagoclearrateobs">
      <value value="1000000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AuRec50primmin">
      <value value="16"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="SaNRec50secomax">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="SaNRxm_mean">
      <value value="0.18"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AuNRgammasecomax">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="EcNRgammaprimmax">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="SaRgammasecomin">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="season">
      <value value="&quot;Continuous variable&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="SaNRgammasecomax">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="EcNRgammasecomin">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AuRgammaprimmax">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="damagerepairrateobs">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="immunegrowthintersectobs">
      <value value="50000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AuNRec50secomin">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="EcRec50primmax">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="maxRinEc">
      <value value="0.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="lungimmunesteadyobs">
      <value value="1000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="EcRxm_mean">
      <value value="0.18"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AuNRxm_sd">
      <value value="0.03"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="EcNRec50secomax">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="SaRgammaprimmax">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="excreterateobs">
      <value value="0.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AuRec50secomax">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="immuneinfluxrateobs">
      <value value="10000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="SaNRec50secomin">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="minRinSa">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AuNRgammasecomin">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AMUstance">
      <value value="&quot;User choice&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AuNRec50primmax">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="SaRec50primmax">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AuRgammaprimmin">
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="SaRxm_mean">
      <value value="0.18"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="SaRec50secomin">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="EcRec50primmin">
      <value value="16"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AuRxm_sd">
      <value value="0.02"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="EcRgammaprimmax">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="SaNRec50primmax">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="EcRgammasecomin">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AuNRgammaprimmax">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="SaRgammaprimmin">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="intakerateobs">
      <value value="0.001"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="SaNRgammaprimmax">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="EcNRgammaprimmin">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="SaNRgammasecomin">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cutoffvalobs">
      <value value="5000000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Ecphagorateobs">
      <value value="1.0E-7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="EcRxm_maxfitcost%">
      <value value="30"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AuNRec50primmin">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="damageincreaserateobs">
      <value value="1.0E-11"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="propSa">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="EcNRec50primmax">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="SaRxm_sd">
      <value value="0.02"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="SaRxm_maxfitcost%">
      <value value="30"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stagephagorateadjustobs">
      <value value="0.7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AuRec50primmax">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="SaNRxm_sd">
      <value value="0.02"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="EcNRec50secomin">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="SaNRec50primmin">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AuNRgammaprimmin">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AuRec50secomin">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="EcNRxm_sd">
      <value value="0.02"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="EcNRgammasecomax">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="EcNRxm_mean">
      <value value="0.18"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AuRgammasecomax">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="removprobobs">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="EcRxm_sd">
      <value value="0.02"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="contamscaleobs">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="patchgrowthpropobs">
      <value value="0.002"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="SaRec50primmin">
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="maxRinSa">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="EcRec50secomax">
      <value value="32"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="EcRgammaprimmin">
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="immunedelayobs">
      <value value="6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="SaRgammasecomax">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AuRxm_maxfitcost%">
      <value value="30"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="phagoinhibobs">
      <value value="1.0E-13"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="surrounddelayobs">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="minRinEc">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Saphagorateobs">
      <value value="5.0E-5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="SaNRgammaprimmin">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="contamweightobs">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="lungmechclearrateobs">
      <value value="1000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AuNRec50secomax">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="SaRec50secomax">
      <value value="128"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="startcap">
      <value value="250000000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AuNRxm_mean">
      <value value="0.18"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AuRgammasecomin">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="propEc">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AuRxm_mean">
      <value value="0.18"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="EcRec50secomin">
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="surroundweightobs">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="moveprob">
      <value value="0.3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="EcNRec50primmin">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="EcRgammasecomax">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="RinAu">
      <value value="0.2"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="final_strategies_experiment5.5" repetitions="200" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="2000"/>
    <metric>deadcounter</metric>
    <metric>sum [mass] of chickens</metric>
    <metric>sum prim-tracker</metric>
    <metric>sum seco-tracker</metric>
    <metric>sum (sublist prim-tracker (ticks - withdrawmoment) (ticks))</metric>
    <metric>sum (sublist seco-tracker (ticks - withdrawmoment) (ticks))</metric>
    <enumeratedValueSet variable="AMUstance">
      <value value="&quot;User choice&quot;"/>
      <value value="&quot;No intervention&quot;"/>
      <value value="&quot;Responsive quick stop&quot;"/>
      <value value="&quot;Responsive complete&quot;"/>
      <value value="&quot;Prophylactic&quot;"/>
      <value value="&quot;Clean responsive&quot;"/>
      <value value="&quot;Regular clean&quot;"/>
      <value value="&quot;Responsive quick stop + clean&quot;"/>
      <value value="&quot;Responsive complete + reg. clean&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="surroundscaleobs">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="lungphagoclearrateobs">
      <value value="1000000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AuRec50primmin">
      <value value="16"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="SaNRec50secomax">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="SaNRxm_mean">
      <value value="0.18"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AuNRgammasecomax">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="SaRgammasecomin">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="EcNRgammaprimmax">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="SaNRgammasecomax">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="season">
      <value value="&quot;Continuous variable&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="EcNRgammasecomin">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="damagerepairrateobs">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AuRgammaprimmax">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="immunegrowthintersectobs">
      <value value="50000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AuNRec50secomin">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="EcRec50primmax">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="maxRinEc">
      <value value="0.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="lungimmunesteadyobs">
      <value value="1000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="EcRxm_mean">
      <value value="0.18"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AuNRxm_sd">
      <value value="0.03"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="EcNRec50secomax">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="SaRgammaprimmax">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AuRec50secomax">
      <value value="50"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="excreterateobs">
      <value value="0.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="immuneinfluxrateobs">
      <value value="10000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="SaNRec50secomin">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="minRinSa">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AuNRgammasecomin">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="SaRec50primmax">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AuNRec50primmax">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AuRgammaprimmin">
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="SaRxm_mean">
      <value value="0.18"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="SaRec50secomin">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="EcRec50primmin">
      <value value="16"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AuRxm_sd">
      <value value="0.02"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="EcRgammaprimmax">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="SaNRec50primmax">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="EcRgammasecomin">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AuNRgammaprimmax">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="SaRgammaprimmin">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="intakerateobs">
      <value value="0.001"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="SaNRgammaprimmax">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="EcNRgammaprimmin">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="SaNRgammasecomin">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cutoffvalobs">
      <value value="5000000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Ecphagorateobs">
      <value value="1.0E-7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="EcRxm_maxfitcost%">
      <value value="30"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AuNRec50primmin">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="damageincreaserateobs">
      <value value="1.0E-11"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="propSa">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="EcNRec50primmax">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="SaRxm_sd">
      <value value="0.02"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="SaRxm_maxfitcost%">
      <value value="30"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="stagephagorateadjustobs">
      <value value="0.7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AuRec50primmax">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="SaNRxm_sd">
      <value value="0.02"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="EcNRec50secomin">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="SaNRec50primmin">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AuNRgammaprimmin">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AuRec50secomin">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="EcNRxm_sd">
      <value value="0.02"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="lunggrowthpropobs">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="EcNRgammasecomax">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="infectprobsobs">
      <value value="0.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="EcNRxm_mean">
      <value value="0.18"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AuRgammasecomax">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="removprobobs">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="EcRxm_sd">
      <value value="0.02"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="patchgrowthpropobs">
      <value value="0.002"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="contamscaleobs">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="SaRec50primmin">
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="maxRinSa">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="EcRec50secomax">
      <value value="32"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="infectscalesobs">
      <value value="6000000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="EcRgammaprimmin">
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="immunedelayobs">
      <value value="6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AuRxm_maxfitcost%">
      <value value="30"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="SaRgammasecomax">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="phagoinhibobs">
      <value value="1.0E-13"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="surrounddelayobs">
      <value value="12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="minRinEc">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Saphagorateobs">
      <value value="5.0E-5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="SaNRgammaprimmin">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="contamweightobs">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="lungmechclearrateobs">
      <value value="1000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="SaRec50secomax">
      <value value="128"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AuNRec50secomax">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="startcap">
      <value value="250000000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AuNRxm_mean">
      <value value="0.18"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AuRgammasecomin">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="propEc">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="AuRxm_mean">
      <value value="0.18"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="surroundweightobs">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="EcRec50secomin">
      <value value="8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="seasonadjustobs">
      <value value="-20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="moveprob">
      <value value="0.3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="EcNRec50primmin">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="EcRgammasecomax">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="RinAu">
      <value value="0.2"/>
    </enumeratedValueSet>
  </experiment>
</experiments>
@#$#@#$#@
@#$#@#$#@
default
0.0
-0.2 0 0.0 1.0
0.0 1 1.0 0.0
0.2 0 0.0 1.0
link direction
true
0
Line -7500403 true 150 150 90 180
Line -7500403 true 150 150 210 180
@#$#@#$#@
0
@#$#@#$#@
