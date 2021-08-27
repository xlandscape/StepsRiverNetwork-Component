"""
Script for documenting the code of the StepsRiverNetwork component.
"""
import os
import base.documentation
import StepsRiverNetwork

root_folder = os.path.abspath(os.path.join(os.path.dirname(base.__file__), ".."))
base.documentation.document_component(
    StepsRiverNetwork.StepsRiverNetwork("StepsRiverNetwork", None, None),
    os.path.join(root_folder, "..", "variant", "StepsRiverNetwork", "README.md"),
    os.path.join(root_folder, "..", "variant", "mc.xml")
)
base.documentation.write_changelog(
    "StepsRiverNetwork component",
    StepsRiverNetwork.StepsRiverNetwork.VERSION,
    os.path.join(root_folder, "..", "variant", "StepsRiverNetwork", "CHANGELOG.md")
)
base.documentation.write_contribution_notes(
    os.path.join(root_folder, "..", "variant", "StepsRiverNetwork", "CONTRIBUTING.md"))
