
Project overview - ${project.project_name} (${project.customer_project_reference})
----------------------------------------------------------------------------

:Project name: ${project.project_name}
:Project id: ${project.project_id}
:Customer reference: ${project.customer_project_reference}
:Order received: ${project.order_received()}
:Contract received: ${project.contract_received()}
:Samples received: ${project.samples_received()}
:Queue date: ${project.queued_date()}
:Application: ${project.application}
:Best practice bioinformatice: ${project.best_practice_bioinformatics}
:Sequencing lanes ordered: ${project.sequencing_units_ordered()}
:Number of samples: ${project.no_of_samples}
:UPPNEX project id: ${project.uppnex_id}
:Project path: /proj/${project.uppnex_id}/INBOX/${project.project_name}

${project.project_sample_name_table()}

${project.project_sample_status_table()}

${project.project_flowcell_summary_table()}


