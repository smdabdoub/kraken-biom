FROM python:3.6.15-buster

WORKDIR /data

RUN pip install biom-format==2.1.10

COPY . .

RUN python setup.py install

CMD ["bash"]
