TAG = r_docker:4.2.1

img:
	docker build --tag=$(TAG) .

shell:
	docker run --rm -it $(TAG) bash

save:
	docker save $(TAG) | gzip > $(shell echo $(TAG) | sed "s/:/-/").tar.gz
