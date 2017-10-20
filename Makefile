ECR_REGISTRY_URL=$(AWS_ACCOUNT_ID).dkr.ecr.$(AWS_DEFAULT_REGION).amazonaws.com

build:
	docker build --build-arg GOOGLE_CLIENT_ID=$GOOGLE_CLIENT_ID --build-arg GOOGLE_CLIENT_SECRET=$GOOGLE_CLIENT_SECRET --build-arg SECRET_KEY=$SECRET_KEY -t cellxgene .

deploy:
	eval $$(aws ecr get-login --region $(AWS_DEFAULT_REGION))
	docker push "$(ECR_REGISTRY_URL)/cellxgene:latest"
