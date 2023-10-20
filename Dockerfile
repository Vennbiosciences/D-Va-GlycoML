FROM nvidia/cuda:11.7.1-base-ubuntu22.04

# install Python 3.10 and pip
RUN : \
    && apt-get update \
    && apt-get install -y \
         python3.10-dev \
         python3-pip \
    && ln -s /usr/bin/python3.10 /usr/local/bin/python \
    && pip install --no-cache --upgrade \
         pip \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# install poetry
RUN pip install poetry

# install poetry dependencies
WORKDIR /app
ADD poetry.lock* pyproject.toml ./
RUN poetry install --no-root

ADD run_train_test.sh ./
ADD ./src ./src

ENTRYPOINT ["/app/run_train_test.sh"]
