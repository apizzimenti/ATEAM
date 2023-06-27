
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

        started = datetime.fromtimestamp(self.start)
        finished = datetime.fromtimestamp(self.end)
        timeToCompletion = timedelta(milliseconds=self.end-self.start)

        self.metadata["started"] = started.strftime("%Y-%m-%d %H:%M:%S")
        self.metadata["completed"] = finished.strftime("%Y-%m-%d %H:%M:%S")
        self.metadata["timeToCompletion"] = str(timeToCompletion)
        self.metadata["steps"] = str(self.chain.steps)
        
        with open("./metadata.json", "w") as w: json.dump(self.metadata, w, indent=4)
