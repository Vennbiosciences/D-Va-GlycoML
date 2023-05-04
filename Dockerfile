FROM artifactory.sfo.venn.bio/docker/library/python:3.8-slim-buster

RUN python -m pip install pipenv

WORKDIR /app
ADD Pipfile* ./

RUN pipenv install --dev --system --deploy

ADD run_train_test.sh ./
ADD ./src ./src

ENTRYPOINT ["/app/run_train_test.sh"]
