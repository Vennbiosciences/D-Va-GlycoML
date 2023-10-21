ARTIFACT ?= intervenn-d-va
TAG ?= latest
IMAGE ?= $(ARTIFACT):$(TAG)
NVIDIA ?= $(shell docker info | grep -q 'Runtimes:.*nvidia' && echo "--runtime=nvidia --gpus=all")

# Folders that should exist in the ./data folder
TRAIN ?= Testing-03-N-HCD-30/
TEST ?= Testing-03-N-HCD-30/
MODEL ?= Models/test-model/


.PHONY: docker
docker:
	docker build \
		-t $(IMAGE) \
		.

.PHONY: run
run:
	sudo rm -rf data/$(MODEL)
	docker run --rm -it $(NVIDIA) \
		--volume $(PWD)/data/:/data/ \
		$(IMAGE) \
		/data/$(TRAIN) \
		/data/$(TEST) \
		/data/$(MODEL)

.PHONY: push
push:
	docker push $(IMAGE)
