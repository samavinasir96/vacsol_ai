"""
URL configuration for vacsol_ai project.

The `urlpatterns` list routes URLs to views. For more information please see:
    https://docs.djangoproject.com/en/4.2/topics/http/urls/
Examples:
Function views
    1. Add an import:  from my_app import views
    2. Add a URL to urlpatterns:  path('', views.home, name='home')
Class-based views
    1. Add an import:  from other_app.views import Home
    2. Add a URL to urlpatterns:  path('', Home.as_view(), name='home')
Including another URLconf
    1. Import the include() function: from django.urls import include, path
    2. Add a URL to urlpatterns:  path('blog/', include('blog.urls'))
"""
from django.contrib import admin
from django.urls import include, path
from .views import *
from .import views
from django.conf.urls.static import static
from django.conf import settings


urlpatterns = [
    path('admin/', admin.site.urls),
    path('',views.home,name='index'),
    path('upload/',views.upload_sequence,name='upload'),
    path('results/',views.get_results,name='results'),
    path('calculate/', views.calculate_features, name='calculate_features'),
    path('get_latest_logs/', views.get_latest_logs, name='get_latest_logs'),
    path('download_csv/', views.download_csv, name='download_csv'),
    path('delete_files/', views.delete_files, name='delete_files'),
    path('features_annotate/', views.features_annotate, name='features_annotate'),
    path('predict_results/', views.predict_results, name='predict_results'),
    path('sequence_upload/', views.sequence_upload, name='sequence_upload'),
    path('delete_files_vacsolai/', views.delete_files_vacsolai, name='delete_files_vacsolai'),

    #path('tutorials/', views.tutorials, name='tutorials'),
    #path('faqs/', views.faqs, name='faqs'),
    #path('contact_us/', views.contact_us, name='contact_us'),
    #path('login_register/', views.login_register, name='login_register'),
] + static(settings.STATIC_URL, document_root=settings.STATIC_ROOT)