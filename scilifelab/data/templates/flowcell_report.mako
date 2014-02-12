
Flowcell summary - ${flowcell.Flowcell}
---------------------------------------

:Date: ${flowcell.Date}
:Flowcell: ${flowcell.Flowcell}
:Position: ${flowcell.FCPosition}
:Instrument: ${flowcell.ScannerID}
:Run mode: ${flowcell.RunMode}
:Setup: ${flowcell.run_setup}
:Index setup: ${flowcell.index_setup()}

${flowcell.flowcell_lane_summary_table(project.flowcell_lanes.get(flowcell.name))}

${flowcell.flowcell_lane_sample_summary_table(project.flowcell_lanes.get(flowcell.name),[s.scilife_name for s in project.project_samples])}

