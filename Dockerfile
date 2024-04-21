FROM python:3.12-slim

WORKDIR /app

COPY . .

RUN mkdir /static

RUN pip install -r requirements.txt
RUN python manage.py collectstatic --noinput

EXPOSE 8000

CMD ["gunicorn", "orbital_calculator.wsgi", "-b", "0.0.0.0:8000"]
