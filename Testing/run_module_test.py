"""Headless driver: runs OMEZarrTest inside Slicer and exits with its status."""

import traceback

import slicer


def main():
    try:
        import OMEZarr

        test = OMEZarr.OMEZarrTest()
        test.runTest()
        print("OMEZARR_TEST_RESULT: PASS")
        status = 0
    except Exception:
        traceback.print_exc()
        print("OMEZARR_TEST_RESULT: FAIL")
        status = 1
    slicer.util.exit(status)


main()
