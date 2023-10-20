ARTIFACT ?= docker/venn/d-va
TAG ?= latest
IMAGE ?= $(ARTIFACT):$(TAG)

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
	rm -rf data/$(MODEL)
	docker run --rm -it \
		--runtime=nvidia \
		--gpus all \
		--volume $(PWD)/data/:/data/ \
		$(IMAGE) \
		/data/$(TRAIN) \
		/data/$(TEST) \
		/data/$(MODEL)

.PHONY: push
push:
	docker push $(IMAGE)
