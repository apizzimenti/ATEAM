
import time
import json
from datetime import datetime, timedelta
from pathlib import Path


class Metadata:
    def __init__(self, chain, directory="./output/statistics/"):
        # Get the experiment metadata.
        with open("./metadata.json") as f: metadata = json.loads(f.read())
        self.metadata = metadata
        self.metadata["field"] = chain.lattice.field.order
        self.metadata["lattice"] = chain.lattice.corners

        self.out = Path(directory)/""
        self.chain = chain

    def __enter__(self):
        self.start = time.time()
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        self.end = time.time()

        # Get some useful runtime statistics.
        diff = self.end-self.start
        diffSeconds = diff/1000
        started = datetime.fromtimestamp(self.start)
        finished = datetime.fromtimestamp(self.end)
        timeToCompletion = timedelta(milliseconds=diff)

        # Compute some stuff; if the rate is over 1, then we're doing more steps
        # than seconds; otherwise, we're doing more seconds than steps.
        rate = self.chain.steps * (1/diffSeconds)
        if rate >= 1: self.metadata["stepsPerSecond"] = rate
        else: self.metadata["secondsPerStep"] = 1/rate

        self.metadata["started"] = started.strftime("%Y-%m-%d %H:%M:%S")
        self.metadata["completed"] = finished.strftime("%Y-%m-%d %H:%M:%S")
        self.metadata["timeToCompletion"] = str(timeToCompletion)
        self.metadata["steps"] = str(self.chain.steps)
        
        with open("./metadata.json", "w") as w: json.dump(self.metadata, w, indent=4)
