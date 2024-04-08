
o2physics_add_dpl_workflow(pid1
                  SOURCES pid1.cxx
                  PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                  COMPONENT_NAME Practice)

o2physics_add_dpl_workflow(r2p2-4-id
                  SOURCES r2p2-4-id.cxx
                  PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                  COMPONENT_NAME Practice)

o2physics_add_dpl_workflow(cros-r2p2-4-id
                  SOURCES cros-r2p2-4-id.cxx
                  PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                  COMPONENT_NAME Practice)
o2physics_add_dpl_workflow(mcpid
                  SOURCES MCPID.cxx
                  PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                  COMPONENT_NAME Practice)
o2physics_add_dpl_workflow(test
                  SOURCES test.cxx
                  PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                  COMPONENT_NAME Practice)
o2physics_add_dpl_workflow(slfr2
                  SOURCES slfr2.cxx
                  PUBLIC_LINK_LIBRARIES O2Physics::AnalysisCore
                  COMPONENT_NAME Practice)


