import subprocess
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from lovis4u.DataProcessing import BedGraph
from lovis4u.Manager import lovis4uError


class BigWigCommandTests(unittest.TestCase):
    def test_bigwig_paths_are_passed_as_separate_process_arguments(self):
        with tempfile.TemporaryDirectory(prefix="lovis4u synthetic output ") as output_dir:
            output_path = str(Path(output_dir) / "synthetic coverage.bedgraph")
            tool_path = "/synthetic tool path/bigWigToBedGraph"
            input_path = "/synthetic input path/coverage profile.bw"

            def write_synthetic_coverage(command, check):
                self.assertEqual(command, [tool_path, input_path, output_path])
                self.assertTrue(check)
                Path(output_path).write_text(
                    "synthetic_contig\t0\t4\t2.5\n", encoding="utf-8"
                )

            parameters = SimpleNamespace(
                args={"bigWigToBedGraph_binary": tool_path}
            )
            fake_temp_file = SimpleNamespace(name=output_path)
            with patch(
                "lovis4u.DataProcessing.tempfile.NamedTemporaryFile",
                return_value=fake_temp_file,
            ):
                with patch(
                    "lovis4u.DataProcessing.subprocess.run",
                    side_effect=write_synthetic_coverage,
                ) as run:
                    coverage = BedGraph(
                        input_path,
                        {"synthetic_contig": 4},
                        parameters,
                        isbigWig=True,
                    )

            run.assert_called_once_with(
                [tool_path, input_path, output_path], check=True
            )
            self.assertEqual(
                coverage.coverage["synthetic_contig"].tolist(),
                [2.5, 2.5, 2.5, 2.5],
            )

    def test_conversion_failure_is_reported(self):
        with tempfile.TemporaryDirectory(prefix="lovis4u synthetic output ") as output_dir:
            output_path = str(Path(output_dir) / "synthetic coverage.bedgraph")
            parameters = SimpleNamespace(
                args={"bigWigToBedGraph_binary": "/synthetic tool/bigWigToBedGraph"}
            )
            fake_temp_file = SimpleNamespace(name=output_path)
            failure = subprocess.CalledProcessError(1, ["synthetic-tool"])

            with patch(
                "lovis4u.DataProcessing.tempfile.NamedTemporaryFile",
                return_value=fake_temp_file,
            ):
                with patch(
                    "lovis4u.DataProcessing.subprocess.run",
                    side_effect=failure,
                ) as run:
                    with self.assertRaises(lovis4uError):
                        BedGraph(
                            "/synthetic input/coverage.bw",
                            {"synthetic_contig": 4},
                            parameters,
                            isbigWig=True,
                        )

            run.assert_called_once_with(
                [
                    "/synthetic tool/bigWigToBedGraph",
                    "/synthetic input/coverage.bw",
                    output_path,
                ],
                check=True,
            )


if __name__ == "__main__":
    unittest.main()
